import sys, os, base64, json, traceback, logging, asyncio, hashlib
from pathlib import Path
from functools import partial
from datetime import datetime, timezone
from fastapi import FastAPI, UploadFile, File, Form, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements
from tos.forensic_artifacts.pdf_generator import generate_v21_dossier
from tos.governance.station_sop import (
    LAW_CANON, LAW_CANON_HASH, STATION_METADATA, TOTAL_LAWS,
    DETERMINISTIC_COUNT, HEURISTIC_COUNT,
    LAW_METHOD_CLASSIFICATIONS, SCORE_DEFINITIONS,
    ARCHITECTURE_WEIGHTS, BAYESIAN_FORMULA
)
from tos.generation.dispatcher import GenerationDispatcher
from tos.enrichment.gemini_compiler import get_compiler
from tos.utils.type_guards import force_bytes
from tos.nkg.manager import get_nkg
from tos.governance.modality_matrix import resolve_method, compute_matrix_hash, get_matrix_meta
from tos.engine.adjudicator import adjudicate_laws, AdjudicationInput
from tos.schemas.response_v1 import ToscaniniResponse
from tos.engine.batch_processor import process_batch
from tos.schemas.batch_v1 import BatchResponse, BatchStructureResult, BatchSummary, FailingLawCount
from tos.security.auth import verify_api_key, validate_upload, validate_batch_upload
from tos.telemetry.usage_logger import log_usage as log_usage_telemetry

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


def _compute_coord_hash(structure) -> str:
    """PIL-NOT-18: Coordinate-Only Hashing (Institutional Standard)."""
    coord_string = "".join(
        f"{round(a.pos[0],3)}{round(a.pos[1],3)}{round(a.pos[2],3)}"
        for a in structure.atoms
    )
    return hashlib.sha256(coord_string.encode()).hexdigest()


def _run_physics_sync(content_bytes: bytes, candidate_id: str, mode: str, t3_category: str):
    # 🛡️ PIL-IMM-20: Immutable Ingestion
    structure = IngestionProcessor.run(force_bytes(content_bytes), "origin.pdb", "Audit", mode)
    
    # 🛡️ PIL-NOT-18: Coordinate-Only Hashing (Institutional Standard)
    coord_hash = _compute_coord_hash(structure)

    full_t1, coverage, _ = Tier1Measurements.run_full_audit(structure, user_intent=t3_category)



    # PIL-CAL-02: Detect experimental source
    _is_exp_source = getattr(structure.confidence, "is_experimental", False) if hasattr(structure, 'confidence') else False
    _method_type = getattr(structure.confidence, "method", "predicted") if _is_exp_source else "predicted"

    # Stages 3+4: Adjudication (delegated to pure function)
    adj_result = adjudicate_laws(AdjudicationInput(
        raw_measurements=full_t1,
        coverage=coverage,
        is_experimental=_is_exp_source,
        method_type=_method_type,
        architecture=t3_category,
    ))

    # Unpack result for payload assembly
    res_t1 = adj_result.law_rows
    failing_det = adj_result.failing_deterministic
    verdict = adj_result.verdict
    det_score = adj_result.deterministic_score
    adv_score = adj_result.advisory_score
    det_passed = adj_result.det_passed
    heur_passed = adj_result.heur_passed
    sm = adj_result.strategic_math
    s6_final, w_arch, m_s8 = sm.s6, sm.w_arch, sm.m_s8
    p_score, suppression = sm.p_score, sm.suppression

    conf_source, conf_val, conf_prov = Tier1Measurements.detect_confidence_source(structure)

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
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "governance_fingerprint": {
                "canon_hash": LAW_CANON_HASH,
                "matrix_hash": compute_matrix_hash(),
                "matrix_schema_version": get_matrix_meta()["schema_version"],
                "policy_ref": get_matrix_meta()["policy_ref"],
            },
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
    try: log_usage_telemetry(get_nkg(), payload, "/ingest")
    except: pass
    return payload

@app.post("/ingest", response_model=ToscaniniResponse)
async def ingest(mode: str = Form(...), candidate_id: str = Form(...), t3_category: str = Form("NONE"), file: UploadFile = File(None), _auth=Depends(verify_api_key)):
    if file: content = await file.read()
    else: content, _, _ = GenerationDispatcher.acquire(candidate_id, None)
    return await asyncio.get_event_loop().run_in_executor(None, partial(_run_physics_sync, content, candidate_id, mode, t3_category))



@app.post("/v1/batch", response_model=BatchResponse)
async def batch_ingest(
    mode: str = Form("Audit"),
    t3_category: str = Form("NONE"),
    file: UploadFile = File(...),
    _auth=Depends(verify_api_key),
):
    """
    Process a ZIP of PDB files and return per-structure verdicts with summary.
    Synchronous execution — results returned directly.
    """
    zip_bytes = await file.read()

    try:
        batch_result = await asyncio.get_event_loop().run_in_executor(
            None,
            partial(process_batch, zip_bytes, _run_physics_sync, mode, t3_category),
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    # Convert dataclass results to Pydantic response
    response_results = []
    for r in batch_result.results:
        response_results.append(BatchStructureResult(
            filename=r.filename,
            candidate_id=r.candidate_id,
            success=r.success,
            response=r.payload if r.success else None,
            error=r.error,
        ))

    summary = BatchSummary(
        total=batch_result.summary.total,
        passed=batch_result.summary.passed,
        vetoed=batch_result.summary.vetoed,
        indeterminate=batch_result.summary.indeterminate,
        errors=batch_result.summary.errors,
        mean_deterministic_score=batch_result.summary.mean_deterministic_score,
        common_failing_laws=[
            FailingLawCount(**fl) for fl in batch_result.summary.common_failing_laws
        ],
    )

    # Log usage per structure
    try:
        nkg = get_nkg()
        for r in batch_result.results:
            if r.success and r.payload:
                log_usage_telemetry(nkg, r.payload, "/v1/batch")
    except: pass

    return BatchResponse(results=response_results, summary=summary)

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
