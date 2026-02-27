"""
Toscanini Phase B2 — Celery Tasks (Week 7: State-Tracked)
"""
import os
import logging
import requests
from celery import Celery
from celery.utils.log import get_task_logger
from worker.job_state import (
    update_job,
    STATE_RUNNING, STATE_SUCCESS, STATE_FAILED
)

logger          = get_task_logger(__name__)
REDIS_URL       = os.environ.get("REDIS_URL", "redis://redis:6379/0")
BRAIN_URL       = os.environ.get("BRAIN_URL", "http://brain:8000")
MAX_JOB_SECONDS = 1800

app = Celery("toscanini_worker", broker=REDIS_URL, backend=REDIS_URL)
app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
    task_track_started=True,
    task_time_limit=MAX_JOB_SECONDS,
    task_soft_time_limit=MAX_JOB_SECONDS - 60,
    worker_prefetch_multiplier=1,
)


@app.task(bind=True, name="worker.tasks.execute_openmm", max_retries=1)
def execute_openmm(self, job_id: str, pdb_bytes_hex: str, protocol_config: dict):
    """Execute OpenMM equilibration on GPU."""
    update_job(job_id, state=STATE_RUNNING, logs="Starting OpenMM...")
    logger.info(f"[{job_id}] OpenMM starting")

    try:
        import openmm as mm
        import openmm.app as app_mm
        import openmm.unit as unit
        from pdbfixer import PDBFixer
        import io

        pdb_bytes     = bytes.fromhex(pdb_bytes_hex)
        temperature_k = protocol_config.get("temperature_k", 300)
        sim_steps     = protocol_config.get("sim_steps", 1000000)
        timestep_fs   = protocol_config.get("timestep_fs", 2)

        update_job(job_id, logs="Fixing structure and adding hydrogens...")
        fixer = PDBFixer(pdbfile=io.StringIO(pdb_bytes.decode("utf-8")))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        update_job(job_id, logs="Building system...")
        forcefield = app_mm.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        modeller   = app_mm.Modeller(fixer.topology, fixer.positions)
        modeller.addSolvent(forcefield, model="tip3p", padding=10 * unit.angstroms)

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app_mm.PME,
            nonbondedCutoff=1.0 * unit.nanometers,
            constraints=app_mm.HBonds,
        )
        integrator = mm.LangevinMiddleIntegrator(
            temperature_k * unit.kelvin,
            1.0 / unit.picosecond,
            timestep_fs * unit.femtoseconds,
        )

        try:
            platform = mm.Platform.getPlatformByName("CUDA")
            update_job(job_id, logs="Using CUDA GPU")
        except Exception:
            platform = mm.Platform.getPlatformByName("CPU")
            update_job(job_id, logs="Using CPU (no GPU)")

        simulation = app_mm.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)

        update_job(job_id, logs="Energy minimization...")
        simulation.minimizeEnergy(maxIterations=1000)

        update_job(job_id, logs=f"Running {sim_steps} steps...")
        simulation.step(sim_steps)

        positions = simulation.context.getState(getPositions=True).getPositions()
        out = io.StringIO()
        app_mm.PDBFile.writeFile(simulation.topology, positions, out)
        refined_pdb = out.getvalue()

        update_job(job_id, state=STATE_SUCCESS,
                   logs=f"OpenMM complete: {sim_steps} steps")

        return {
            "status":          "success",
            "job_id":          job_id,
            "protocol":        "openmm",
            "refined_pdb_hex": refined_pdb.encode("utf-8").hex(),
            "logs":            f"OpenMM complete: {sim_steps} steps"
        }

    except Exception as e:
        update_job(job_id, state=STATE_FAILED, error=str(e))
        logger.error(f"[{job_id}] OpenMM failed: {e}")
        return {"status": "failed", "job_id": job_id, "error": str(e)}


@app.task(bind=True, name="worker.tasks.execute_rosetta", max_retries=1)
def execute_rosetta(self, job_id: str, pdb_bytes_hex: str, xml_protocol: str):
    """Execute Rosetta FastRelax."""
    update_job(job_id, state=STATE_RUNNING, logs="Starting Rosetta...")
    logger.info(f"[{job_id}] Rosetta starting")

    ROSETTA_BIN = os.environ.get(
        "ROSETTA_BIN",
        "/opt/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease"
    )

    if not os.path.exists(ROSETTA_BIN):
        msg = "Rosetta binary not found"
        update_job(job_id, state=STATE_FAILED, error=msg)
        return {"status": "failed", "job_id": job_id, "error": msg}

    try:
        import subprocess, tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "input.pdb")
            xml_path = os.path.join(tmpdir, "protocol.xml")
            out_path = os.path.join(tmpdir, "output.pdb")

            with open(pdb_path, "wb") as f:
                f.write(bytes.fromhex(pdb_bytes_hex))
            with open(xml_path, "w") as f:
                f.write(xml_protocol)

            update_job(job_id, logs="Running Rosetta FastRelax...")
            result = subprocess.run(
                f"{ROSETTA_BIN} -in:file:s {pdb_path} "
                f"-parser:protocol {xml_path} -out:file:o {out_path} "
                f"-out:nstruct 1 -out:overwrite -ex1 -ex2 "
                f"-use_input_sc -mute all -unmute protocols.relax",
                shell=True, capture_output=True, text=True,
                timeout=MAX_JOB_SECONDS - 120, cwd=tmpdir
            )

            if result.returncode != 0:
                raise RuntimeError(f"Rosetta failed: {result.stderr[:300]}")
            if not os.path.exists(out_path):
                raise RuntimeError("Rosetta produced no output PDB")

            with open(out_path, "rb") as f:
                refined_pdb = f.read()

            update_job(job_id, state=STATE_SUCCESS, logs="Rosetta complete")
            return {
                "status":          "success",
                "job_id":          job_id,
                "protocol":        "rosetta",
                "refined_pdb_hex": refined_pdb.hex(),
                "logs":            result.stdout[-500:]
            }

    except Exception as e:
        update_job(job_id, state=STATE_FAILED, error=str(e))
        return {"status": "failed", "job_id": job_id, "error": str(e)}


@app.task(name="worker.tasks.auto_callback")
def auto_callback(execution_result: dict, original_audit_id: str,
                  user_email: str = None):
    """Post-execution: send refined PDB to brain for re-audit."""
    job_id = execution_result.get("job_id", "UNKNOWN")

    if execution_result.get("status") != "success":
        logger.error(f"[{job_id}] Execution failed — skipping callback")
        return {"status": "skipped", "reason": "execution_failed"}

    try:
        response = requests.post(
            f"{BRAIN_URL}/refinement/internal-callback",
            json={
                "original_audit_id": original_audit_id,
                "refined_pdb_hex":   execution_result["refined_pdb_hex"],
                "job_id":            job_id,
                "protocol":          execution_result.get("protocol"),
                "user_email":        user_email
            },
            timeout=120
        )

        if response.status_code == 200:
            res             = response.json()
            refined_id      = res.get("refined_audit_id")
            comparison_url  = res.get("comparison_url")
            update_job(job_id,
                       refined_audit_id=refined_id,
                       comparison_url=comparison_url)
            logger.info(f"[{job_id}] Callback success → {refined_id}")
            return {
                "status":           "success",
                "job_id":           job_id,
                "refined_audit_id": refined_id,
                "comparison_url":   comparison_url
            }
        else:
            err = response.text[:200]
            update_job(job_id, state=STATE_FAILED, error=err)
            return {"status": "failed", "error": err}

    except Exception as e:
        update_job(job_id, state=STATE_FAILED, error=str(e))
        return {"status": "failed", "error": str(e)}
