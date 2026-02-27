"""
Toscanini Phase B2 — Celery Task Definitions
GPU worker tasks for managed refinement execution.

Tasks:
- execute_openmm: OpenMM equilibration
- execute_rosetta: Rosetta FastRelax
- auto_callback: Post-execution re-audit trigger
"""
import os
import time
import tempfile
import logging
import requests
from celery import Celery
from celery.utils.log import get_task_logger

logger = get_task_logger(__name__)

# Redis connection
REDIS_URL = os.environ.get("REDIS_URL", "redis://redis:6379/0")
BRAIN_URL = os.environ.get("BRAIN_URL", "http://brain:8000")
MAX_JOB_SECONDS = 1800  # 30 min hard timeout

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
    worker_prefetch_multiplier=1,  # One job at a time (GPU resource)
)


@app.task(bind=True, name="worker.tasks.execute_openmm", max_retries=1)
def execute_openmm(self, job_id: str, pdb_bytes_hex: str, protocol_config: dict):
    """
    Execute OpenMM equilibration on GPU.

    Args:
        job_id: Unique job identifier
        pdb_bytes_hex: PDB file contents as hex string
        protocol_config: Protocol parameters (steps, temperature, etc.)

    Returns:
        dict with status, refined_pdb_hex, logs
    """
    logger.info(f"[{job_id}] Starting OpenMM execution")

    try:
        import openmm as mm
        import openmm.app as app_mm
        import openmm.unit as unit
        from pdbfixer import PDBFixer
        import io

        pdb_bytes = bytes.fromhex(pdb_bytes_hex)

        # Protocol parameters
        temperature_k = protocol_config.get("temperature_k", 300)
        sim_steps = protocol_config.get("sim_steps", 1000000)  # 2ns default
        timestep_fs = protocol_config.get("timestep_fs", 2)

        logger.info(f"[{job_id}] Protocol: {sim_steps} steps @ {temperature_k}K")

        # Fix structure
        fixer = PDBFixer(pdbfile=io.StringIO(pdb_bytes.decode("utf-8")))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        # Create system
        forcefield = app_mm.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        modeller = app_mm.Modeller(fixer.topology, fixer.positions)
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

        # Use GPU if available
        try:
            platform = mm.Platform.getPlatformByName("CUDA")
            logger.info(f"[{job_id}] Using CUDA platform")
        except Exception:
            platform = mm.Platform.getPlatformByName("CPU")
            logger.info(f"[{job_id}] Falling back to CPU platform")

        simulation = app_mm.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)

        # Minimize
        logger.info(f"[{job_id}] Running energy minimization")
        simulation.minimizeEnergy(maxIterations=1000)

        # Run equilibration
        logger.info(f"[{job_id}] Running {sim_steps} steps equilibration")
        simulation.step(sim_steps)

        # Extract refined structure
        positions = simulation.context.getState(getPositions=True).getPositions()
        output = io.StringIO()
        app_mm.PDBFile.writeFile(simulation.topology, positions, output)
        refined_pdb = output.getvalue()

        logger.info(f"[{job_id}] OpenMM execution complete")

        return {
            "status": "success",
            "job_id": job_id,
            "protocol": "openmm",
            "refined_pdb_hex": refined_pdb.encode("utf-8").hex(),
            "steps_completed": sim_steps,
            "logs": f"OpenMM equilibration complete: {sim_steps} steps"
        }

    except Exception as e:
        logger.error(f"[{job_id}] OpenMM execution failed: {str(e)}")
        return {
            "status": "failed",
            "job_id": job_id,
            "protocol": "openmm",
            "error": str(e),
            "logs": f"OpenMM failed: {str(e)}"
        }


@app.task(bind=True, name="worker.tasks.execute_rosetta", max_retries=1)
def execute_rosetta(self, job_id: str, pdb_bytes_hex: str, xml_protocol: str):
    """
    Execute Rosetta FastRelax.

    Args:
        job_id: Unique job identifier
        pdb_bytes_hex: PDB file contents as hex string
        xml_protocol: Rosetta XML protocol string

    Returns:
        dict with status, refined_pdb_hex, logs
    """
    import subprocess
    import tempfile

    logger.info(f"[{job_id}] Starting Rosetta execution")

    ROSETTA_BIN = os.environ.get(
        "ROSETTA_BIN",
        "/opt/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease"
    )

    if not os.path.exists(ROSETTA_BIN):
        return {
            "status": "failed",
            "job_id": job_id,
            "protocol": "rosetta",
            "error": "Rosetta binary not found. Check ROSETTA_BIN env var.",
            "logs": "Rosetta not installed"
        }

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_bytes = bytes.fromhex(pdb_bytes_hex)

            # Write input files
            pdb_path = os.path.join(tmpdir, "input.pdb")
            xml_path = os.path.join(tmpdir, "protocol.xml")
            out_path = os.path.join(tmpdir, "output.pdb")

            with open(pdb_path, "wb") as f:
                f.write(pdb_bytes)
            with open(xml_path, "w") as f:
                f.write(xml_protocol)

            # Build Rosetta command
            cmd = [
                ROSETTA_BIN,
                f"-in:file:s {pdb_path}",
                f"-parser:protocol {xml_path}",
                f"-out:file:o {out_path}",
                "-out:nstruct 1",
                "-out:overwrite",
                "-ex1", "-ex2",
                "-use_input_sc",
                "-mute all",
                "-unmute protocols.relax"
            ]

            logger.info(f"[{job_id}] Running Rosetta FastRelax")
            result = subprocess.run(
                " ".join(cmd),
                shell=True,
                capture_output=True,
                text=True,
                timeout=MAX_JOB_SECONDS - 120,
                cwd=tmpdir
            )

            if result.returncode != 0:
                raise RuntimeError(f"Rosetta failed: {result.stderr[:500]}")

            # Read output
            if not os.path.exists(out_path):
                raise RuntimeError("Rosetta produced no output PDB")

            with open(out_path, "rb") as f:
                refined_pdb = f.read()

            logger.info(f"[{job_id}] Rosetta execution complete")

            return {
                "status": "success",
                "job_id": job_id,
                "protocol": "rosetta",
                "refined_pdb_hex": refined_pdb.hex(),
                "logs": result.stdout[-1000:]
            }

    except Exception as e:
        logger.error(f"[{job_id}] Rosetta execution failed: {str(e)}")
        return {
            "status": "failed",
            "job_id": job_id,
            "protocol": "rosetta",
            "error": str(e),
            "logs": f"Rosetta failed: {str(e)}"
        }


@app.task(name="worker.tasks.auto_callback")
def auto_callback(execution_result: dict, original_audit_id: str, user_email: str = None):
    """
    Post-execution: upload refined structure to brain for re-audit.
    Triggered automatically after execute_openmm or execute_rosetta.

    Args:
        execution_result: Result dict from execution task
        original_audit_id: Original audit ID for comparison
        user_email: Optional email for notification
    """
    job_id = execution_result.get("job_id", "UNKNOWN")

    if execution_result.get("status") != "success":
        logger.error(f"[{job_id}] Execution failed — skipping callback")
        return {"status": "skipped", "reason": "execution_failed"}

    logger.info(f"[{job_id}] Triggering auto-callback to brain")

    try:
        refined_pdb_hex = execution_result["refined_pdb_hex"]
        refined_pdb_bytes = bytes.fromhex(refined_pdb_hex)

        # Call brain internal callback endpoint
        response = requests.post(
            f"{BRAIN_URL}/refinement/internal-callback",
            json={
                "original_audit_id": original_audit_id,
                "refined_pdb_hex": refined_pdb_hex,
                "job_id": job_id,
                "protocol": execution_result.get("protocol"),
                "user_email": user_email
            },
            timeout=120
        )

        if response.status_code == 200:
            result = response.json()
            logger.info(f"[{job_id}] Auto-callback success: {result.get('refined_audit_id')}")
            return {
                "status": "success",
                "job_id": job_id,
                "refined_audit_id": result.get("refined_audit_id"),
                "comparison_url": result.get("comparison_url")
            }
        else:
            logger.error(f"[{job_id}] Auto-callback failed: {response.text[:200]}")
            return {"status": "failed", "error": response.text[:200]}

    except Exception as e:
        logger.error(f"[{job_id}] Auto-callback error: {str(e)}")
        return {"status": "failed", "error": str(e)}
