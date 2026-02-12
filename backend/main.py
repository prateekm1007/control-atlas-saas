import sys, os, base64, json, traceback, logging
from pathlib import Path
from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements
from tos.engine.tier3_predict import StrategicPredictor
from tos.saas_engine.selector import ArchitectureSelector
from tos.forensic_artifacts.pdf_generator import generate_v21_dossier
from tos.governance.station_sop import LAW_CANON, SCORE_DEFINITIONS
from tos.glossary.law_glossary import list_all_law_ids, get_law_explanation
from tos.discovery.resolver import DiscoveryResolver
from tos.generation.dispatcher import GenerationDispatcher
from tos.enrichment.gemini_compiler import get_compiler
from tos.nkg.manager import get_nkg
from tos.governance.determinism import unify_architecture_authority
from tos.utils.type_guards import force_bytes
import tos.glossary.epistemic_definitions as epistemic_mod

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("toscanini.brain")

app = FastAPI(title="Toscanini OS", version="21.35.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

@app.post("/ingest")
async def ingest(mode: str = Form(...), candidate_id: str = Form(...), t3_category: str = Form("NONE"), file: UploadFile = File(None)):
    try:
        content = await file.read() if file else GenerationDispatcher.acquire(candidate_id, None)[0]
        structure = IngestionProcessor.run(force_bytes(content), "origin.pdb", "Audit", mode)
        
        full_t1 = Tier1Measurements.run_full_audit(structure)
        res_t1 = []
        for lid in list_all_law_ids():
            res_t1.append({
                "law_id": lid, 
                "status": full_t1[lid][0], 
                "measurement": str(full_t1[lid][1]), 
                "title": LAW_CANON[lid]['title'], 
                "principle": LAW_CANON[lid]['principle']
            })
        
        phys_score = int((sum(1 for l in res_t1 if l['status'] == "PASS") / 11) * 100)
        actual_conf = round(structure.confidence.mean_plddt, 1)
        t3_res = StrategicPredictor.calculate(phys_score, actual_conf, 0.85)
        arch_res = unify_architecture_authority(ArchitectureSelector.select(structure, user_intent=t3_category))

        ctx = {"v": "PASS" if phys_score==100 else "VETO", "s": phys_score, "p": t3_res['probability'], "c": actual_conf, "arch": arch_res['authoritative_category']}
        narratives = get_compiler().synthesize_dossier_content(ctx)

        payload = {
            "verdict": {"binary": ctx["v"], "physical_score": phys_score, "confidence_score": actual_conf},
            "definitions": epistemic_mod.DEFINITIONS,
            "governance": {"program_id": f"PRG-{structure.audit_id[:8]}", "state": "SEALED"},
            "provenance": {"source": candidate_id, "audit_id": structure.audit_id, "hash": structure.get_coordinate_hash()},
            "tier1": {"laws": res_t1}, "tier3": t3_res, "architecture": arch_res,
            "witness_reports": narratives,
            "pdb_b64": base64.b64encode(force_bytes(content)).decode(),
            "visualization": {"ca_coords": [list(a.pos) for a in structure.atoms if a.atom_name == 'CA']}
        }
        get_nkg().record_audit(payload)
        pdf_b = generate_v21_dossier(payload)
        return {**payload, "pdf_b64": base64.b64encode(pdf_b).decode()}
    except Exception as e:
        logger.error(traceback.format_exc())
        return JSONResponse(status_code=500, content={"verdict": {"binary": "ERROR"}, "error": str(e)})

@app.get("/definitions")
def get_definitions(): return epistemic_mod.DEFINITIONS
@app.post("/search")
async def search(query: str = Form(...)): return DiscoveryResolver.resolve(query)
@app.get("/nkg/records")
def read_nkg(limit: int = 20): return get_nkg().read_records(limit)
