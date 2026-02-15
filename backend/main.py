import sys, os, base64, json, traceback, logging, asyncio, hashlib
from pathlib import Path
from functools import partial
from datetime import datetime, timezone
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
    LAW_METHOD_CLASSIFICATIONS, SCORE_DEFINITIONS,
    ARCHITECTURE_WEIGHTS, BAYESIAN_FORMULA
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
    import numpy as _np
    if isinstance(obj, dict): return {k: _sanitize_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)): return [_sanitize_for_json(v) for v in obj]
    elif isinstance(obj, (_np.integer,)): return int(obj)
    elif isinstance(obj, (_np.floating,)): return float(obj)
    elif isinstance(obj, (_np.bool_,)): return bool(obj)
    elif isinstance(obj, _np.ndarray): return obj.tolist()
    return obj

def _compute_verdict(res_t1, coverage, fatal_fringe):
    """🛡️ PIL-LCY-19: 3-State Forward-Only State Machine"""
    cov_threshold = float(LAW_CANON["LAW-105"]["threshold"])
    if coverage < cov_threshold: return "INDETERMINATE"
    det_failures = [r for r in res_t1 if r["method"] == "deterministic" and r["status"] in ("FAIL", "VETO")]
    if det_failures: return "VETO"
    if fatal_fringe: return "INDETERMINATE"
    heur_failures = [r for r in res_t1 if r["method"] == "heuristic" and r["status"] in ("FAIL", "VETO")]
    if heur_failures: return "INDETERMINATE"
    return "PASS"

def _run_physics_sync(content_bytes: bytes, candidate_id: str, mode: str, t3_category: str):
    # 🛡️ PIL-IMM-20: Immutable Ingestion
    structure = IngestionProcessor.run(force_bytes(content_bytes), "origin.pdb", "Audit", mode)
    
    # 🛡️ PIL-NOT-18: Coordinate-Only Hashing (Institutional Standard)
    coord_string = "".join([f"{round(a.pos[0],3)}{round(a.pos[1],3)}{round(a.pos[2],3)}" for a in structure.atoms])
    coord_hash = hashlib.sha256(coord_string.encode()).hexdigest()

    full_t1, coverage, fatal_fringe = Tier1Measurements.run_full_audit(structure, user_intent=t3_category)

    res_t1, failing_det = [], []
    for lid in sorted(LAW_CANON.keys()):
        m = full_t1.get(lid, {"observed": 0, "deviation": "0.0", "sample": 0, "status": "FAIL"})
        row = {
            "law_id": lid, "title": LAW_CANON[lid]["title"], "status": m["status"],
            "method": LAW_METHOD_CLASSIFICATIONS.get(lid, "unknown"),
            "observed": m["observed"], "threshold": LAW_CANON[lid]["threshold"],
            "operator": LAW_CANON[lid]["operator"], "units": LAW_CANON[lid]["unit"],
            "deviation": m["deviation"], "sample_size": m["sample"],
            "scope": LAW_CANON[lid]["scope"], "principle": LAW_CANON[lid].get("principle", "N/A")
        }
        res_t1.append(row)
        if m["status"] in ("FAIL", "VETO") and row["method"] == "deterministic":
            failing_det.append(f"{lid}: {LAW_CANON[lid]['title']}")

    det_passed = sum(1 for r in res_t1 if r["status"] == "PASS" and r["method"] == "deterministic")
    det_score = int((det_passed / max(DETERMINISTIC_COUNT, 1)) * 100)
    heur_passed = sum(1 for r in res_t1 if r["status"] == "PASS" and r["method"] == "heuristic")
    adv_score = int((heur_passed / max(HEURISTIC_COUNT, 1)) * 100)

    conf_source, conf_val, conf_prov = Tier1Measurements.detect_confidence_source(structure)
    verdict = _compute_verdict(res_t1, coverage, fatal_fringe)

    # 🛡️ PIL-MTH-12: Strategic Success Math
    raw_s6 = round(det_score / 100.0, 4)
    w_arch = ARCHITECTURE_WEIGHTS.get(t3_category, 1.0)
    m_s8 = 0.0 # NKG Penalty
    
    # Coverage Gate Logic (Review Item 5: Explicitly zero index and integrity on failure)
    cov_threshold = float(LAW_CANON["LAW-105"]["threshold"])
    if coverage < cov_threshold:
        s6_final = 0.0
        p_score = 0
        suppression = "Coverage gate enforces zeroing of prioritization index."
    else:
        s6_final = raw_s6
        p_score = int(s6_final * w_arch * (1.0 - m_s8) * 100)
        suppression = None

    payload = _sanitize_for_json({
        "verdict": {
            "binary": verdict, "deterministic_score": det_score, "advisory_score": adv_score,
            "physical_score": det_score, "confidence_score": conf_val, 
            "det_passed": det_passed, "det_total": DETERMINISTIC_COUNT, 
            "heur_passed": heur_passed, "heur_total": HEURISTIC_COUNT, 
            "coverage_pct": round(coverage, 2), "suppression_reason": suppression
        },
        "governance": {
            "audit_id": str(hashlib.md5(coord_hash.encode()).hexdigest()[:8]).upper(),
            "station_version": STATION_METADATA["version"],
            "timestamp_utc": datetime.now(timezone.utc).isoformat()
        },
        "provenance": {"source": candidate_id, "hash": coord_hash, "byte_count": len(content_bytes)},
        "tier1": {"laws": res_t1}, "tier3": {"probability": p_score},
        "characterization": Tier1Measurements.compute_structural_characterization(structure),
        "confidence_meta": {
            "source_type": conf_source, "provenance_method": conf_prov,
            **Tier1Measurements.decompose_confidence(structure)
        },
        "strategic_math": {
            "s6": s6_final, "w_arch": w_arch, "m_s8": m_s8, "architecture": t3_category
        },
        "witness_reports": get_compiler().synthesize_dossier_content({
            "v": verdict, "s": det_score, "c": coverage, "arch": t3_category, "killer_laws": failing_det
        }),
        "ai_model_used": get_compiler().model_used,
        "pdb_b64": base64.b64encode(content_bytes).decode(),
    })

    try:
        payload["pdf_b64"] = base64.b64encode(generate_v21_dossier(payload)).decode()
    except Exception:
        logger.error(f"PDF FAILED: {traceback.format_exc()}")
        payload["pdf_b64"] = ""

    try: get_nkg().record_audit(payload)
    except: pass
    return payload

@app.post("/ingest")
async def ingest(mode: str = Form(...), candidate_id: str = Form(...), t3_category: str = Form("NONE"), file: UploadFile = File(None)):
    if file: content = await file.read()
    else: content, _, _ = GenerationDispatcher.acquire(candidate_id, None)
    return await asyncio.get_event_loop().run_in_executor(None, partial(_run_physics_sync, content, candidate_id, mode, t3_category))

@app.get("/health")
def health(): return {"status": "operational", "version": STATION_METADATA["version"]}

@app.get("/laws")
def get_laws(): return {"laws": LAW_CANON, "total_laws": TOTAL_LAWS, "deterministic_count": DETERMINISTIC_COUNT}

@app.get("/nkg")
def get_nkg_records():
    try:
        nkg = get_nkg()
        records = nkg.read_records(limit=50)
        return {"vetoes": records.get("vetoes", []), "successes": records.get("successes", []), 
                "total_vetoes": len(records.get("vetoes", [])), "total_successes": len(records.get("successes", []))}
    except: return {"vetoes": [], "successes": []}

@app.get("/definitions")
def get_definitions(): return SCORE_DEFINITIONS

@app.post("/search")
async def search(query: str = Form(...)):
    from tos.discovery.resolver import DiscoveryResolver
    results = DiscoveryResolver.resolve(query)
    return {"results": results, "count": len(results)}
