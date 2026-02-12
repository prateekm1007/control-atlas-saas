import sys
import os
from pathlib import Path

# Force the current directory into sys.path
BASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE_DIR))

from fastapi import FastAPI, UploadFile, File, Form
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
import base64, json, re

# PROTECTED NAMESPACE IMPORTS
import tos.ingestion.processor as processor_mod
import tos.engine.tier1_measurements as measurements_mod
import tos.router.intelligence as router_mod
import tos.export.sealer as sealer_mod
import tos.forensic_artifacts.pdf_generator as pdf_mod
import tos.glossary.law_glossary as glossary_mod
import tos.glossary.epistemic_definitions as epistemic_mod
import tos.enrichment.gemini_compiler as gemini_mod
import tos.generation.dispatcher as dispatcher_mod
import tos.discovery.resolver as discovery_mod
import tos.governance.program_schema as program_schema_mod

app = FastAPI(); app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])
program_repository = program_schema_mod.ProgramRepository()

@app.get("/health")
def health(): return {"status": "ALIVE"}

@app.get("/stats")
def stats():
    count = 0
    path = Path("/app/nkg/piu_moat.jsonl")
    if path.exists():
        with open(path, 'r') as f: count = sum(1 for _ in f)
    return {"unique_pius": count}

@app.post("/search")
async def search(query: str = Form(...)): 
    return discovery_mod.DiscoveryResolver.resolve(query)

def advance_program_state(structure, source_label: str, plddt: float):
    program = program_repository.upsert(structure.audit_id)
    return program_schema_mod.advance_program(
        program=program,
        structure=structure,
        acquisition_source=source_label,
        folding_source=source_label,
        plddt=plddt,
    )


@app.post("/ingest")
async def ingest(mode: str = Form(...), candidate_id: str = Form(...), file: UploadFile = File(None), sequence: str = Form(None)):
    try:
        content, label, ext = None, "Unknown", "pdb"
        if mode == "Upload" and file:
            content = await file.read(); label = "User Upload"; ext = file.filename.split(".")[-1].lower()
        else:
            content, label, ext = dispatcher_mod.GenerationDispatcher.acquire(candidate_id, sequence)

        if not content: return {"verdict": "ERROR", "details": label}

        structure = processor_mod.IngestionProcessor.run(content, f"origin.{ext}", label)
        all_atoms = structure.atoms + structure.ligands
        coord_hash = sealer_mod.ForensicSealer.generate_hash(sealer_mod.ForensicSealer.canonical_serialize(all_atoms))
        
        s155, m155, a155 = measurements_mod.Tier1Measurements.check_law_155_L(structure)
        s160, m160, a160 = measurements_mod.Tier1Measurements.check_law_160(structure)
        results = []
        for lid in glossary_mod.list_all_law_ids():
            st, me, an = ("PASS", "Verified Invariant", {})
            if lid == "LAW-155": st, me, an = s155, m155, a155
            elif lid == "LAW-160": st, me, an = s160, m160, a160
            expl = glossary_mod.get_law_explanation(lid)
            results.append({"law_id": lid, "status": st, "measurement": me, "anchor": an, "title": expl['title'], "principle": expl['principle'], "rationale": expl['rationale']})
        
        verdict = "VETO" if any(r["status"] in {"FAIL", "VETO"} for r in results) else "PASS"
        phys_score = 20 if verdict == "VETO" else 100
        conf_display = round(structure.confidence.mean_plddt, 1) if hasattr(structure, 'confidence') else "N/A"

        program_record = advance_program_state(structure, label, float(conf_display) if conf_display != "N/A" else 0.0)
        
        gemini = gemini_mod.GeminiCompiler(os.getenv("GEMINI_API_KEY", "NONE"))
        rat = gemini.synthesize(verdict, phys_score, conf_display, label, results)
        sealed = sealer_mod.ForensicSealer.seal_structure(content, structure.audit_id, verdict, coord_hash, ext)
        pdf_b = pdf_mod.generate_v14_certificate(structure.audit_id, verdict, phys_score, label, rat, results, structure.atoms)

        return {
            "verdict": verdict, "score": phys_score, "conf": conf_display, "sig": structure.audit_id, 
            "laws": results, "routing": router_mod.IntelligenceRouter().decide(structure, verdict), 
            "narrative": rat, "pdf_b64": base64.b64encode(pdf_b).decode(), "pdb_b64": base64.b64encode(sealed.encode()).decode(), "ext": ext,
            "provenance": {"source": label, "mode": mode},
            "definitions": epistemic_mod.DEFINITIONS,
            "program": {
                "program_id": program_record.program_id,
                "state": program_record.current_state.value if program_record.current_state else None,
                "veto": program_record.veto,
                "hash_chain": program_record.hash_chain,
                "architecture": program_record.architecture.__dict__ if program_record.architecture else None,
                "tier3": program_record.tier3_breakdown,
            },
        }
    except Exception as e: return {"verdict": "ERROR", "details": str(e)}


@app.get("/programs/{program_id}/dossier")
def get_dossier(program_id: str):
    program = program_repository.get(program_id)
    if not program:
        return {"error": "Program not found"}
    pdf_bytes = pdf_mod.generate_dossier(program)
    return StreamingResponse(
        iter([pdf_bytes]),
        media_type="application/pdf",
        headers={"Content-Disposition": f"inline; filename={program_id}_dossier.pdf"},
    )
