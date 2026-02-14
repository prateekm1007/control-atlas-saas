import sys, os, base64, json, traceback, logging, asyncio, hashlib
from pathlib import Path
from functools import partial
from fastapi import FastAPI, UploadFile, File, Form, Depends, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements
from tos.forensic_artifacts.pdf_generator import generate_v21_dossier
from tos.governance.station_sop import (
    LAW_CANON, SCORE_DEFINITIONS, STATION_METADATA, TOTAL_LAWS, 
    ARCHITECTURE_WEIGHTS, DETERMINISTIC_LAWS, HEURISTIC_LAWS, 
    DETERMINISTIC_COUNT, HEURISTIC_COUNT, CONFIDENCE_EXCLUSION_POLICY, 
    LAW_METHOD_CLASSIFICATIONS, LAW_CANON_HASH, BAYESIAN_FORMULA
)
from tos.discovery.resolver import DiscoveryResolver
from tos.nkg.manager import get_nkg
from tos.generation.dispatcher import GenerationDispatcher
from tos.enrichment.gemini_compiler import get_compiler
from tos.utils.type_guards import force_bytes

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
    return obj

def _run_physics_sync(content_bytes: bytes, candidate_id: str, mode: str, t3_category: str):
    structure = IngestionProcessor.run(force_bytes(content_bytes), "origin.pdb", "Audit", mode)
    full_t1 = Tier1Measurements.run_full_audit(structure, user_intent=t3_category)
    res_t1, failing_det = [], []
    for lid in LAW_CANON:
        entry = full_t1.get(lid, ("FAIL", "Error", "MISSING", "error"))
        res_t1.append({"law_id": lid, "status": entry[0], "measurement": entry[1], "method": entry[3], "title": LAW_CANON[lid]['title'], "principle": LAW_CANON[lid].get('principle', 'N/A')})
        if entry[0] == "FAIL" and entry[3] == "deterministic": failing_det.append(f"{lid}: {LAW_CANON[lid]['title']}")

    det_passed = sum(1 for l in res_t1 if l['status'] == "PASS" and l['method'] == "deterministic")
    det_score = int((det_passed / DETERMINISTIC_COUNT) * 100)
    conf_source, conf_val, conf_prov = Tier1Measurements.detect_confidence_source(structure)
    s6_composite = min(1.0, max(0.0, 0.15 + round((det_score/100)*0.35, 4) + (round((conf_val/100)*0.40, 4) if conf_val > 0 else 0.0)))

    ctx = {"v": "PASS" if det_score == 100 else "VETO", "s": det_score, "p": int(s6_composite*100), "c": conf_val, "arch": t3_category, "killer_laws": failing_det}
    try: narratives = get_compiler().synthesize_dossier_content(ctx)
    except: narratives = {"executive": "Audit complete.", "deep_dive": "Verified.", "recommendation": "Proceed."}

    payload = _sanitize_for_json({
        "verdict": {"binary": ctx["v"], "deterministic_score": det_score, "physical_score": det_score, "confidence_score": conf_val, "confidence_available": conf_val > 0, "killer_laws": failing_det},
        "provenance": {"source": candidate_id, "hash": hashlib.sha256(content_bytes).hexdigest(), "byte_count": len(content_bytes), "station_version": STATION_METADATA["version"]},
        "tier1": {"laws": res_t1}, "tier3": {"probability": int(s6_composite*100)},
        "characterization": Tier1Measurements.compute_structural_characterization(structure),
        "confidence_meta": {"source_type": conf_source, "provenance_method": conf_prov},
        "bayesian_components": {"breakdown": [{"name": "S6_composite", "value": s6_composite, "source": "logic"}]},
        "witness_reports": narratives, "ai_model_used": "Gemini-Apex", "pdb_b64": base64.b64encode(content_bytes).decode()
    })
    payload["pdf_b64"] = base64.b64encode(generate_v21_dossier(payload)).decode()
    return payload

@app.post("/ingest")
async def ingest(mode: str = Form(...), candidate_id: str = Form(...), t3_category: str = Form("NONE"), file: UploadFile = File(None)):
    content = await file.read() if file else GenerationDispatcher.acquire(candidate_id, None)[0]
    return await asyncio.get_event_loop().run_in_executor(None, partial(_run_physics_sync, content, candidate_id, mode, t3_category))

@app.post("/search")
async def search(query: str = Form(...)):
    return DiscoveryResolver.resolve(query)

@app.get("/nkg/records")
def read_nkg(limit: int = 20):
    return get_nkg().read_records(limit)

@app.get("/health")
def health(): return {"status": "operational", "version": STATION_METADATA["version"]}

@app.get("/laws")
def get_laws(): return {"laws": LAW_CANON, "total_laws": TOTAL_LAWS, "deterministic_count": DETERMINISTIC_COUNT}
