import sys, os, base64, json, traceback, logging, asyncio, hashlib
from pathlib import Path
from functools import partial
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements
from tos.forensic_artifacts.pdf_generator import generate_v21_dossier
from tos.governance.station_sop import (
    LAW_CANON, STATION_METADATA, TOTAL_LAWS,
    DETERMINISTIC_COUNT, HEURISTIC_COUNT,
    LAW_METHOD_CLASSIFICATIONS
)
from tos.generation.dispatcher import GenerationDispatcher
from tos.enrichment.gemini_compiler import get_compiler
from tos.utils.type_guards import force_bytes
from tos.nkg.manager import get_nkg

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("toscanini.brain")

app = FastAPI(title="Toscanini OS", version=STATION_METADATA["version"])
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


def _sanitize_for_json(obj):
    """Recursively convert numpy types to native Python for JSON."""
    import numpy as _np
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_sanitize_for_json(v) for v in obj]
    elif isinstance(obj, (_np.integer,)):
        return int(obj)
    elif isinstance(obj, (_np.floating,)):
        return float(obj)
    elif isinstance(obj, (_np.bool_,)):
        return bool(obj)
    elif isinstance(obj, _np.ndarray):
        return obj.tolist()
    return obj


def _compute_verdict(res_t1, coverage, fatal_fringe):
    """
    3-state verdict:
      PASS:          All deterministic laws pass, coverage >= 30%
      VETO:          Any deterministic core law fails
      INDETERMINATE: Coverage < 30%, or only fringe/heuristic failures
    """
    # Gate 1: Coverage insufficient — cannot trust results
    if coverage < 30.0:
        return "INDETERMINATE"

    # Gate 2: Any deterministic core failure → hard VETO
    det_failures = [
        r for r in res_t1
        if r["method"] == "deterministic"
        and r["status"] in ("FAIL", "VETO")
    ]
    if det_failures:
        return "VETO"

    # Gate 3: Fatal fringe escalation → soft INDETERMINATE
    if fatal_fringe:
        return "INDETERMINATE"

    # Gate 4: Heuristic-only failures → INDETERMINATE (not VETO)
    heur_failures = [
        r for r in res_t1
        if r["method"] == "heuristic"
        and r["status"] in ("FAIL", "VETO")
    ]
    if heur_failures:
        return "INDETERMINATE"

    return "PASS"


def _run_physics_sync(content_bytes: bytes, candidate_id: str, mode: str, t3_category: str):
    """Core physics pipeline — returns complete audit payload."""
    structure = IngestionProcessor.run(force_bytes(content_bytes), "origin.pdb", "Audit", mode)

    # Engine returns 3-tuple: (results_dict, coverage_pct, fatal_fringe_bool)
    full_t1, coverage, fatal_fringe = Tier1Measurements.run_full_audit(
        structure, user_intent=t3_category
    )

    # Assemble law results list
    res_t1 = []
    failing_det = []
    for lid in sorted(LAW_CANON.keys()):
        entry = full_t1.get(lid, ("FAIL", "Missing implementation", "MISSING", "error"))
        status = entry[0]
        method = entry[3] if len(entry) > 3 else LAW_METHOD_CLASSIFICATIONS.get(lid, "unknown")
        row = {
            "law_id": lid,
            "status": status,
            "measurement": str(entry[1]),
            "method": method,
            "title": LAW_CANON[lid]["title"],
            "principle": LAW_CANON[lid].get("principle", "N/A"),
        }
        res_t1.append(row)
        if status in ("FAIL", "VETO") and method == "deterministic":
            failing_det.append(f"{lid}: {LAW_CANON[lid]['title']}")

    # Scores
    det_passed = sum(1 for r in res_t1 if r["status"] == "PASS" and r["method"] == "deterministic")
    det_score = int((det_passed / max(DETERMINISTIC_COUNT, 1)) * 100)
    heur_passed = sum(1 for r in res_t1 if r["status"] == "PASS" and r["method"] == "heuristic")
    adv_score = int((heur_passed / max(HEURISTIC_COUNT, 1)) * 100)

    # Confidence
    conf_source, conf_val, conf_prov = Tier1Measurements.detect_confidence_source(structure)

    # 3-state verdict
    verdict = _compute_verdict(res_t1, coverage, fatal_fringe)

    # Simplified composite for pilot
    s6_composite = (det_score / 100.0) * 0.85
    if conf_val > 0:
        s6_composite += (conf_val / 100.0) * 0.15

    # AI narratives — best-effort, never crash the audit
    try:
        compiler = get_compiler()
        ctx = {
            "v": verdict, "s": det_score, "p": int(s6_composite * 100),
            "c": conf_val, "arch": t3_category,
            "killer_laws": failing_det, "failing_laws": [],
        }
        narratives = compiler.synthesize_dossier_content(ctx)
        ai_model = compiler.model_used
    except Exception:
        narratives = {
            "executive": "Audit complete.",
            "deep_dive": "Verified.",
            "recommendation": "Proceed per metrics.",
        }
        ai_model = "Internal (deterministic fallback)"

    # Assemble payload
    payload = _sanitize_for_json({
        "verdict": {
            "binary": verdict,
            "deterministic_score": det_score,
            "advisory_score": adv_score,
            "physical_score": det_score,
            "confidence_score": conf_val,
            "confidence_available": conf_val > 0,
            "killer_laws": failing_det,
            "laws_passed": det_passed + heur_passed,
            "laws_total": TOTAL_LAWS,
            "det_passed": det_passed,
            "det_total": DETERMINISTIC_COUNT,
            "heur_passed": heur_passed,
            "heur_total": HEURISTIC_COUNT,
            "coverage_pct": round(coverage, 1),
        },
        "provenance": {
            "source": candidate_id,
            "hash": hashlib.sha256(content_bytes).hexdigest(),
            "byte_count": len(content_bytes),
            "station_version": STATION_METADATA["version"],
        },
        "tier1": {"laws": res_t1},
        "tier3": {"probability": int(s6_composite * 100)},
        "characterization": Tier1Measurements.compute_structural_characterization(structure),
        "confidence_meta": {
            "source_type": conf_source,
            "provenance_method": conf_prov,
        },
        "bayesian_components": {
            "breakdown": [
                {"name": "S6_composite", "value": round(s6_composite, 4), "source": "logic"},
            ]
        },
        "witness_reports": narratives,
        "ai_model_used": ai_model,
        "pdb_b64": base64.b64encode(content_bytes).decode(),
    })

    # PDF — never crash the audit over a PDF failure
    try:
        payload["pdf_b64"] = base64.b64encode(generate_v21_dossier(payload)).decode()
    except Exception as e:
        logger.error(f"PDF generation failed: {e}")
        payload["pdf_b64"] = ""

    # NKG: Record every audit (successes AND failures)
    try:
        get_nkg().record_audit(payload)
    except Exception as e:
        logger.warning(f"NKG recording failed: {e}")

    return payload


# ── Routes ──────────────────────────────────────────────────────
@app.post("/ingest")
async def ingest(
    mode: str = Form(...),
    candidate_id: str = Form(...),
    t3_category: str = Form("NONE"),
    file: UploadFile = File(None),
):
    try:
        if file:
            content = await file.read()
            if not content:
                raise HTTPException(status_code=400, detail="Empty file uploaded")
        else:
            content, _, _ = GenerationDispatcher.acquire(candidate_id, None)

        result = await asyncio.get_event_loop().run_in_executor(
            None, partial(_run_physics_sync, content, candidate_id, mode, t3_category)
        )
        return result

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Audit failed: {traceback.format_exc()}")
        return JSONResponse(
            status_code=500,
            content={
                "error": "Audit processing failed",
                "detail": str(e),
                "verdict": {"binary": "ERROR"},
            },
        )


@app.get("/health")
def health():
    return {"status": "operational", "version": STATION_METADATA["version"]}


@app.get("/laws")
def get_laws():
    return {
        "laws": LAW_CANON,
        "total_laws": TOTAL_LAWS,
        "deterministic_count": DETERMINISTIC_COUNT,
        "heuristic_count": HEURISTIC_COUNT,
    }


@app.get("/definitions")
def get_definitions():
    from tos.governance.station_sop import SCORE_DEFINITIONS
    return SCORE_DEFINITIONS


@app.post("/search")
async def search(query: str = Form(...)):
    from tos.discovery.resolver import DiscoveryResolver
    results = DiscoveryResolver.resolve(query)
    return {"results": results, "count": len(results)}


@app.get("/nkg")
def get_nkg_records():
    try:
        nkg = get_nkg()
        records = nkg.read_records(limit=50)
        return {
            "vetoes": records.get("vetoes", []),
            "successes": records.get("successes", []),
            "total_vetoes": len(records.get("vetoes", [])),
            "total_successes": len(records.get("successes", [])),
        }
    except Exception as e:
        return {"vetoes": [], "successes": [], "error": str(e)}
