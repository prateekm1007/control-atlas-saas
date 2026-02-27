"""
Toscanini Phase B2 — Celery Tasks (Week 12: Security-Hardened)

Changes from Week 7:
  - Rosetta: shell=True removed, list args used, os.setsid + killpg for full
    process-group kill on timeout or error. No orphaned GPU processes.
  - SoftTimeLimitExceeded handled cleanly before hard SIGKILL.
  - Celery time limits already set (1800s hard / 1740s soft).
"""
import os
import signal
import logging
import subprocess
import tempfile
import requests
from contextlib import contextmanager
from celery import Celery
from celery.exceptions import SoftTimeLimitExceeded
from celery.utils.log import get_task_logger
from worker.job_state import (
    update_job,
    STATE_RUNNING, STATE_SUCCESS, STATE_FAILED
)

logger          = get_task_logger(__name__)
REDIS_URL       = os.environ.get("REDIS_URL", "redis://redis:6379/0")
BRAIN_URL       = os.environ.get("BRAIN_URL", "http://brain:8000")
MAX_JOB_SECONDS = 1800   # hard wall — matches task_time_limit

app = Celery("toscanini_worker", broker=REDIS_URL, backend=REDIS_URL)
app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
    task_track_started=True,
    task_time_limit=MAX_JOB_SECONDS,           # SIGKILL after 30 min
    task_soft_time_limit=MAX_JOB_SECONDS - 60, # SoftTimeLimitExceeded at 29 min
    worker_prefetch_multiplier=1,
    worker_max_tasks_per_child=2,              # recycle worker every 2 jobs
)


# ── Week 12: Hardened subprocess context manager ──────────────────────────────

@contextmanager
def _hard_timeout_process(cmd: list, timeout_sec: int):
    """
    Launch a subprocess in its own process group.
    On timeout or any exception: SIGTERM → 3 s grace → SIGKILL.
    Guarantees no orphaned GPU processes. Ever.

    Usage:
        with _hard_timeout_process(cmd, timeout_sec=900) as proc:
            stdout, stderr = proc.communicate()
    """
    proc = subprocess.Popen(
        cmd,                           # list — no shell injection possible
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        preexec_fn=os.setsid,          # new process group → killpg kills children
    )
    try:
        yield proc
        proc.wait(timeout=timeout_sec)
    except subprocess.TimeoutExpired:
        _kill_group(proc)
        raise
    except Exception:
        _kill_group(proc)
        raise


def _kill_group(proc: subprocess.Popen) -> None:
    """SIGTERM then SIGKILL the entire process group."""
    try:
        pgid = os.getpgid(proc.pid)
        os.killpg(pgid, signal.SIGTERM)
        try:
            proc.wait(timeout=3)
        except subprocess.TimeoutExpired:
            pass
        os.killpg(pgid, signal.SIGKILL)
    except ProcessLookupError:
        pass   # already dead
    try:
        proc.wait()
    except Exception:
        pass


# ── OpenMM Task ───────────────────────────────────────────────────────────────

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

    except SoftTimeLimitExceeded:
        update_job(job_id, state=STATE_FAILED,
                   logs="Soft time limit reached — job terminated cleanly")
        raise
    except Exception as e:
        update_job(job_id, state=STATE_FAILED, error=str(e))
        logger.error(f"[{job_id}] OpenMM failed: {e}")
        return {"status": "failed", "job_id": job_id, "error": str(e)}


# ── Rosetta Task ──────────────────────────────────────────────────────────────

@app.task(bind=True, name="worker.tasks.execute_rosetta", max_retries=1)
def execute_rosetta(self, job_id: str, pdb_bytes_hex: str, xml_protocol: str):
    """
    Execute Rosetta FastRelax.

    Week 12 hardening:
      - shell=True removed; cmd is a list (no injection surface)
      - _hard_timeout_process: os.setsid + killpg kills entire process group
      - SoftTimeLimitExceeded handled before Celery hard-kills
    """
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
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "input.pdb")
            xml_path = os.path.join(tmpdir, "protocol.xml")
            out_path = os.path.join(tmpdir, "output.pdb")

            with open(pdb_path, "wb") as f:
                f.write(bytes.fromhex(pdb_bytes_hex))
            with open(xml_path, "w") as f:
                f.write(xml_protocol)

            # ── Week 12: list args, no shell=True ─────────────────────────
            cmd = [
                ROSETTA_BIN,
                "-in:file:s",      pdb_path,
                "-parser:protocol", xml_path,
                "-out:file:o",     out_path,
                "-out:nstruct",    "1",
                "-out:overwrite",
                "-ex1",
                "-ex2",
                "-use_input_sc",
                "-mute",           "all",
                "-unmute",         "protocols.relax",
            ]

            update_job(job_id, logs="Running Rosetta FastRelax...")

            rosetta_timeout = MAX_JOB_SECONDS - 120  # 28 min — leave 2 min for callback

            with _hard_timeout_process(cmd, timeout_sec=rosetta_timeout) as proc:
                stdout_bytes, stderr_bytes = proc.communicate(
                    timeout=rosetta_timeout
                )

            stdout_text = stdout_bytes.decode("utf-8", errors="replace")
            stderr_text = stderr_bytes.decode("utf-8", errors="replace")

            if proc.returncode != 0:
                raise RuntimeError(
                    f"Rosetta exited {proc.returncode}. "
                    f"STDERR: {stderr_text[-500:]}"
                )
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
                "logs":            stdout_text[-500:]
            }

    except subprocess.TimeoutExpired:
        update_job(job_id, state=STATE_FAILED,
                   error="Rosetta exceeded hard timeout — process killed")
        return {"status": "failed", "job_id": job_id,
                "error": "Timeout"}

    except SoftTimeLimitExceeded:
        update_job(job_id, state=STATE_FAILED,
                   logs="Soft time limit reached — Rosetta terminated cleanly")
        raise

    except Exception as e:
        update_job(job_id, state=STATE_FAILED, error=str(e))
        return {"status": "failed", "job_id": job_id, "error": str(e)}


# ── Callback Task ─────────────────────────────────────────────────────────────

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
            res            = response.json()
            refined_id     = res.get("refined_audit_id")
            comparison_url = res.get("comparison_url")
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
