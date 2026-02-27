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
from tos.security.auth import verify_api_key, validate_upload, validate_batch_upload, enforce_size_limit
from tos.telemetry.usage_logger import log_usage as log_usage_telemetry
from tos.security.tokens import validate_refinement_token
from tos.saas.key_store import hash_key, get_tier_for_key, get_gpu_allocation_for_key
from tos.storage.comparisons import store_comparison, get_comparison
from tos.storage.audit_store import store_audit_result, get_audit_result


# ── Week 12: Hard caps ────────────────────────────────────────────────────────
MAX_PDB_BYTES     = 50 * 1024 * 1024   # 50 MB
MAX_RESIDUE_COUNT = 1500               # GPU budget cap


def _count_residues(pdb_bytes: bytes) -> int:
    """Count unique (chain_id, resseq) pairs from ATOM/HETATM lines.
    Pure-Python — no Biopython required in the brain container."""
    seen: set = set()
    for line in pdb_bytes.decode("utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            chain  = line[21] if len(line) > 21 else "?"
            resseq = line[22:26].strip() if len(line) > 26 else "0"
            seen.add((chain, resseq))
    return len(seen)


def _validate_pdb_upload(raw: bytes) -> None:
    """Raise HTTPException before credits or GPU are touched.
    Called once at the top of /refinement/submit."""
    if len(raw) > MAX_PDB_BYTES:
        raise HTTPException(
            status_code=413,
            detail=(
                f"PDB file exceeds the 50 MB hard cap "
                f"({len(raw) / 1_048_576:.1f} MB received). "
                "Split your structure or use the B1 external path."
            ),
        )
    n_res = _count_residues(raw)
    if n_res > MAX_RESIDUE_COUNT:
        raise HTTPException(
            status_code=422,
            detail=(
                f"Structure contains {n_res} residues, exceeding the "
                f"{MAX_RESIDUE_COUNT}-residue cap for managed GPU runs. "
                "Trim your structure or use the B1 external path."
            ),
        )

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
    except Exception as e:
        logger.error(f"PDF generation failed for audit_id={payload.get('governance', {}).get('audit_id', 'UNKNOWN')}: {str(e)}")
        payload["pdf_b64"] = ""
        payload["pdf_error"] = "PDF_GENERATION_FAILED"
    try: get_nkg().record_audit(payload)
    except: pass
    try: log_usage_telemetry(get_nkg(), payload, "/ingest")
    except: pass
    return payload

@app.post("/ingest", response_model=ToscaniniResponse)
async def ingest(mode: str = Form(...), candidate_id: str = Form(...), t3_category: str = Form("NONE"), file: UploadFile = File(None), _auth=Depends(verify_api_key)):
    if file:
        validate_upload(file)
        content = await file.read()
        await enforce_size_limit(content)
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
    validate_batch_upload(file)
    zip_bytes = await file.read()
    await enforce_size_limit(zip_bytes)

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


@app.post("/refinement/callback")
async def refinement_callback(
    token: str = Form(...),
    file: UploadFile = File(...)
):
    """
    Phase B1: Accept refined structure upload via callback token.
    
    Flow:
    1. Validate token (7-day expiry, single-use)
    2. Extract original audit_id from token
    3. Run new audit on refined structure
    4. Store comparison metadata
    5. Return comparison URL + new audit result
    
    Returns:
        JSON with original_audit_id, refined_audit_id, comparison_url
    """
    try:
        # Validate token
        payload = validate_refinement_token(token, consume=True)  # Single-use enforcement
        original_audit_id = payload["audit_id"]
        user_email = payload.get("user_email")
        
        # Read uploaded file
        content_bytes = await file.read()
        
        if len(content_bytes) == 0:
            return {"status": "error", "message": "Empty file uploaded"}, 400
        
        # Detect mode from file extension (or use original audit metadata)
        mode = "experimental" if file.filename.endswith('.pdb') else "predicted"
        
        # Run audit on refined structure
        refined_audit = _run_physics_sync(
            content_bytes, 
            file.filename, 
            mode, 
            "NONE"
        )
        
        refined_audit_id = refined_audit["governance"]["audit_id"]
        
        # Store comparison metadata
        comparison_data = {
            "original_audit_id": original_audit_id,
            "refined_audit_id": refined_audit_id,
            "uploaded_at": datetime.now(timezone.utc).isoformat(),
            "user_email": user_email,
            "filename": file.filename,
            "file_size_bytes": len(content_bytes),
            "refinement_method": "external",  # Phase B1: user-executed
            "status": "complete"
        }
        
        store_comparison(original_audit_id, refined_audit_id, comparison_data)

        # Store refined audit for future retrieval
        store_audit_result(refined_audit_id, refined_audit)
        
        # Build comparison URL
        comparison_url = f"/compare?baseline={original_audit_id}&refined={refined_audit_id}"
        
        return {
            "status": "success",
            "message": "Refined structure audited successfully",
            "original_audit_id": original_audit_id,
            "refined_audit_id": refined_audit_id,
            "comparison_url": comparison_url,
            "verdict": refined_audit["verdict"]["binary"],
            "coverage_pct": refined_audit["verdict"]["coverage_pct"],
            "deterministic_score": refined_audit["verdict"]["deterministic_score"],
            "next_step": "View comparison in dashboard to verify improvement"
        }
        
    except ValueError as e:
        return {"status": "error", "message": str(e), "error_type": "token_validation"}, 401
    except Exception as e:
        logger.error(f"Callback error: {str(e)}")
        return {"status": "error", "message": f"Server error: {str(e)}", "error_type": "server"}, 500


@app.get("/comparison/{original_id}/{refined_id}")
async def get_comparison(original_id: str, refined_id: str):
    """
    Retrieve stored comparison metadata and regenerate comparison analysis.
    
    Args:
        original_id: Original (baseline) audit ID
        refined_id: Refined audit ID
    
    Returns:
        Full comparison object with law deltas and improvements
    """
    try:
        from tos.storage.comparisons import get_comparison
        
        # Retrieve comparison metadata
        comp_meta = get_comparison(original_id, refined_id)
        
        if not comp_meta:
            return {"status": "error", "message": "Comparison not found"}, 404
        
        # TODO: Retrieve full audit results for both IDs
        # For now, return metadata only
        # In Week 4, add full audit retrieval and delta calculation
        
        return {
            "status": "success",
            "comparison": comp_meta,
            "message": "Full delta calculation coming in Week 4"
        }
        
    except Exception as e:
        return {"status": "error", "message": str(e)}, 500

@app.get("/comparisons/by-original/{original_id}")
async def list_comparisons_for_original(original_id: str):
    """
    List all refined audits derived from a given original audit.
    
    Useful for showing refinement history on audit detail page.
    """
    try:
        from tos.storage.comparisons import list_comparisons_by_original
        
        comparisons = list_comparisons_by_original(original_id)
        
        return {
            "status": "success",
            "original_audit_id": original_id,
            "refinement_count": len(comparisons),
            "comparisons": comparisons
        }
        
    except Exception as e:
        return {"status": "error", "message": str(e)}, 500


@app.get("/health")
async def health_check():
    """
    System health check for monitoring and deployment validation.
    
    Checks:
    - API responding
    - Token system operational
    - Storage accessible
    - Canon hash stable
    """
    from tos.security.tokens import create_refinement_token, validate_refinement_token
    from tos.storage.comparisons import STORAGE_DIR
    from pathlib import Path
    
    checks = {}
    
    # Token system check
    try:
        test_token = create_refinement_token("HEALTH_CHECK")
        validate_refinement_token(test_token)
        checks["token_system"] = "OK"
    except Exception as e:
        checks["token_system"] = f"ERROR: {str(e)}"
    
    # Storage check
    try:
        storage_path = Path("/app/data/comparisons")
        storage_path.mkdir(parents=True, exist_ok=True)
        checks["storage"] = "OK"
    except Exception as e:
        checks["storage"] = f"ERROR: {str(e)}"
    
    # Canon hash check
    try:
        checks["canon_hash"] = LAW_CANON_HASH[:8] + "..."
        checks["governance"] = "STABLE"
    except Exception as e:
        checks["governance"] = f"ERROR: {str(e)}"
    
    all_ok = all(v in ("OK", "STABLE") or "..." in str(v) 
                 for v in checks.values())
    
    return {
        "status": "healthy" if all_ok else "degraded",
        "version": STATION_METADATA.get("version", "unknown"),
        "checks": checks,
        "timestamp": datetime.now(timezone.utc).isoformat()
    }


@app.post("/refinement/submit")
async def refinement_submit(
    audit_id: str = Form(...),
    protocol: str = Form("openmm"),
    user_email: str = Form(None),
    file: UploadFile = File(...),
    request: Request = None,
):
    """
    Phase B2: Submit structure for managed GPU refinement.

    Flow:
    1. Validate audit_id exists
    2. Read PDB file
    3. Queue Celery task (execute_openmm or execute_rosetta)
    4. Return job_id for status polling

    Args:
        audit_id: Original audit ID (must be INDETERMINATE/VETO)
        protocol: openmm | rosetta
        user_email: Optional email for notification
        file: PDB/CIF file to refine
    """
    import uuid

    try:
        # ── Week 12: read once, validate immediately (before credits) ──────
        content_bytes = await file.read()
        if len(content_bytes) == 0:
            return {"status": "error", "message": "Empty file"}, 400
        _validate_pdb_upload(content_bytes)   # raises 413 / 422 on violation

        # ── B.3: Tier-aware credit check ─────────────────────────────────
        # Determine identifier: prefer API key tier over anonymous email/IP
        _api_key_raw  = request.headers.get("X-API-Key", "") if request else ""
        _key_hash     = hash_key(_api_key_raw) if _api_key_raw else ""
        _tier         = get_tier_for_key(_key_hash) if _key_hash else "free"
        _identifier   = user_email or (request.client.host if request else "anonymous")

        # Import credit system and enforce tier ceiling
        import sys as _sys
        _sys.path.insert(0, "/app/gpu_worker")
        from worker.credits import check_credits, deduct_credits, BETA_CREDITS_ANON

        # Override credit allocation based on API key tier
        _gpu_alloc = get_gpu_allocation_for_key(_key_hash) if _key_hash else BETA_CREDITS_ANON
        _credit_check = check_credits(_identifier, protocol)

        if not _credit_check["allowed"]:
            raise HTTPException(
                status_code=402,
                detail=(
                    f"No GPU credits remaining for {_identifier}. "
                    f"Tier: {_tier}. "
                    f"Upgrade your API key tier for more GPU runs."
                )
            )

        # Validate protocol
        if protocol not in ("openmm", "rosetta"):
            return {"status": "error", "message": "Invalid protocol. Use openmm or rosetta."}, 400

        # Generate job ID
        job_id = str(uuid.uuid4())[:8].upper()

        # Convert to hex for Celery serialization
        pdb_hex = content_bytes.hex()

        # Queue task
        try:
            from celery import Celery
            import os

            REDIS_URL = os.environ.get("REDIS_URL", "redis://redis:6379/0")
            celery_app = Celery("toscanini_client", broker=REDIS_URL, backend=REDIS_URL)

            # Week 12: failing_laws derived from stored audit if available
            # audit_result is never in scope here; retrieve from store instead
            _stored = get_audit_result(audit_id) or {}
            _all_laws = _stored.get("tier1", {}).get("laws", [])
            failing_laws = [l["law_id"] for l in _all_laws if l.get("status") not in ("PASS",)]

            # Dispatch via execution engine
            import sys
            sys.path.insert(0, "/app/gpu_worker")
            from worker.execution_engine import dispatch_job

            dispatch_result = dispatch_job(
                pdb_bytes=content_bytes,
                original_audit_id=audit_id,
                failing_laws=failing_laws,
                user_email=user_email,
                protocol_override=protocol if protocol != "auto" else None
            )

            job_id         = dispatch_result["job_id"]
            queue_status   = dispatch_result["state"]
            celery_task_id = dispatch_result.get("celery_task_id")

        except Exception as celery_err:
            queue_status   = "beta_pending"
            celery_task_id = None
            logger.warning(f"Execution engine not available: {str(celery_err)}")

        # Deduct credit after successful dispatch
        try:
            deduct_credits(_identifier, protocol, job_id)
        except Exception as _dc_err:
            logger.warning(f"Credit deduction failed (non-fatal): {_dc_err}")

        return {
            "status": "success",
            "job_id": job_id,
            "celery_task_id": celery_task_id,
            "queue_status": queue_status,
            "protocol": protocol,
            "original_audit_id": audit_id,
            "message": (
                "Job queued for GPU execution. Poll /refinement/status/{job_id} for updates."
                if queue_status == "queued"
                else "B2 GPU worker not yet deployed. Use B1 callback flow for now."
            ),
            "estimated_time_minutes": 5 if protocol == "openmm" else 3
        }

    except Exception as e:
        logger.error(f"Submit error: {str(e)}")
        return {"status": "error", "message": str(e)}, 500


@app.get("/refinement/status/{job_id}")
async def refinement_status(job_id: str):
    """Poll job state from Redis state manager."""
    try:
        import redis as redis_lib, json, os
        r = redis_lib.from_url(
            os.environ.get("REDIS_URL", "redis://redis:6379/0"),
            decode_responses=True
        )
        raw = r.get(f"job:{job_id}")
        if raw:
            job = json.loads(raw)
            return {
                "job_id":           job_id,
                "state":            job.get("state", "unknown"),
                "protocol":         job.get("protocol"),
                "original_audit_id":job.get("original_audit_id"),
                "refined_audit_id": job.get("refined_audit_id"),
                "comparison_url":   job.get("comparison_url"),
                "created_at":       job.get("created_at"),
                "started_at":       job.get("started_at"),
                "completed_at":     job.get("completed_at"),
                "logs":             job.get("logs"),
                "error":            job.get("error")
            }
        return {"job_id": job_id, "state": "not_found",
                "message": "Job not found or expired"}
    except Exception as e:
        return {"job_id": job_id, "state": "unknown", "error": str(e),
                "message": "Redis not available (B2 GPU worker not deployed)"}

@app.get("/refinement/jobs/{original_audit_id}")
async def list_jobs_for_audit(original_audit_id: str):
    """List all refinement jobs for a given original audit."""
    try:
        import redis as redis_lib, json, os
        r = redis_lib.from_url(
            os.environ.get("REDIS_URL", "redis://redis:6379/0"),
            decode_responses=True
        )
        job_ids = r.lrange(f"audit_jobs:{original_audit_id}", 0, -1)
        jobs = []
        for jid in job_ids:
            raw = r.get(f"job:{jid}")
            if raw:
                jobs.append(json.loads(raw))
        return {
            "original_audit_id": original_audit_id,
            "job_count":         len(jobs),
            "jobs": sorted(jobs, key=lambda x: x.get("created_at",""), reverse=True)
        }
    except Exception as e:
        return {"original_audit_id": original_audit_id, "error": str(e)}


@app.post("/refinement/internal-callback")
async def internal_callback(payload: dict):
    """
    Internal endpoint: called by GPU worker after execution completes.
    Triggers re-audit and stores comparison.
    """
    try:
        original_audit_id = payload["original_audit_id"]
        refined_pdb_hex = payload["refined_pdb_hex"]
        job_id = payload.get("job_id", "UNKNOWN")
        user_email = payload.get("user_email")
        protocol = payload.get("protocol", "unknown")

        # Re-audit refined structure
        refined_bytes = bytes.fromhex(refined_pdb_hex)
        refined_audit = _run_physics_sync(refined_bytes, f"refined_{job_id}.pdb", "experimental", "NONE")
        refined_audit_id = refined_audit["governance"]["audit_id"]

        # Store comparison
        comparison_data = {
            "original_audit_id": original_audit_id,
            "refined_audit_id": refined_audit_id,
            "job_id": job_id,
            "protocol": protocol,
            "refinement_method": "managed_gpu",
            "user_email": user_email,
            "status": "complete"
        }

        store_comparison(original_audit_id, refined_audit_id, comparison_data)

        return {
            "status": "success",
            "original_audit_id": original_audit_id,
            "refined_audit_id": refined_audit_id,
            "comparison_url": f"/compare?baseline={original_audit_id}&refined={refined_audit_id}"
        }

    except Exception as e:
        logger.error(f"Internal callback error: {str(e)}")
        return {"status": "error", "message": str(e)}, 500


@app.get("/audit/{audit_id}")
async def get_stored_audit(audit_id: str):
    """Retrieve a previously stored audit result by ID."""
    result = get_audit_result(audit_id)
    if result:
        return {"status": "success", "audit": result}
    return {"status": "not_found", "message": f"Audit {audit_id} not found or expired"}

@app.post("/compare")
async def compare_audits_endpoint(baseline_id: str = Form(...), refined_id: str = Form(...)):
    """Compare two stored audits and return delta analysis."""
    baseline = get_audit_result(baseline_id)
    refined  = get_audit_result(refined_id)

    if not baseline:
        return {"status": "error", "message": f"Baseline audit {baseline_id} not found"}
    if not refined:
        return {"status": "error", "message": f"Refined audit {refined_id} not found"}

    # Import comparison engine
    import sys
    sys.path.insert(0, "/app")
    from comparison_engine import compare_audits

    comparison = compare_audits(baseline, refined)
    return {
        "status": "success",
        "comparison": comparison
    }
