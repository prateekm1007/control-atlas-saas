"""
Toscanini Phase B2 — Execution Engine
Orchestrates protocol selection and task dispatch.

Wraps Celery task submission with:
- Protocol auto-selection
- Job state initialization
- Error handling
- Estimated time reporting
"""
import os
import uuid
import logging
from typing import Optional, Dict
from worker.protocol_selector import (
    select_protocol, get_openmm_config, get_rosetta_scoreterms,
    PROTOCOL_ROSETTA, PROTOCOL_OPENMM, PROTOCOL_BOTH, PROTOCOL_LOOP
)
from worker.job_state import create_job, STATE_QUEUED

logger    = logging.getLogger("toscanini.execution_engine")
REDIS_URL = os.environ.get("REDIS_URL", "redis://redis:6379/0")


def generate_job_id() -> str:
    """Generate unique 8-char uppercase job ID."""
    return str(uuid.uuid4())[:8].upper()


def dispatch_job(pdb_bytes: bytes,
                 original_audit_id: str,
                 failing_laws: list,
                 user_email: Optional[str] = None,
                 protocol_override: Optional[str] = None) -> Dict:
    """
    Main entry point: select protocol, init job state, dispatch Celery task.

    Args:
        pdb_bytes:          Raw PDB file bytes
        original_audit_id:  Original audit ID (baseline)
        failing_laws:       List of failing law IDs
        user_email:         Optional email for notification
        protocol_override:  Force specific protocol (openmm|rosetta|both)

    Returns:
        dict with job_id, protocol, estimated_minutes, state
    """
    job_id    = generate_job_id()
    pdb_hex   = pdb_bytes.hex()

    # Select protocol
    if protocol_override:
        selection = {
            "protocol":          protocol_override,
            "rationale":         "User-specified protocol",
            "estimated_minutes": 5
        }
    else:
        selection = select_protocol(failing_laws)

    protocol          = selection["protocol"]
    estimated_minutes = selection["estimated_minutes"]

    # Initialize job in Redis
    create_job(job_id, original_audit_id, protocol, user_email)

    logger.info(f"[{job_id}] Dispatching {protocol} job for {original_audit_id}")

    # Dispatch to Celery
    try:
        from celery import Celery
        celery_app = Celery("toscanini_client", broker=REDIS_URL, backend=REDIS_URL)

        if protocol == PROTOCOL_OPENMM:
            config = get_openmm_config(failing_laws)
            task = celery_app.send_task(
                "worker.tasks.execute_openmm",
                args=[job_id, pdb_hex, config],
                queue="refinement",
                link=_make_callback_signature(celery_app, original_audit_id, user_email)
            )

        elif protocol == PROTOCOL_ROSETTA:
            weights  = get_rosetta_scoreterms(failing_laws)
            xml      = _generate_rosetta_xml(original_audit_id, failing_laws, weights)
            task = celery_app.send_task(
                "worker.tasks.execute_rosetta",
                args=[job_id, pdb_hex, xml],
                queue="refinement",
                link=_make_callback_signature(celery_app, original_audit_id, user_email)
            )

        elif protocol in (PROTOCOL_BOTH, PROTOCOL_LOOP):
            # Phase 1: dispatch primary task, secondary chained after callback
            config  = get_openmm_config(failing_laws)
            task = celery_app.send_task(
                "worker.tasks.execute_openmm",
                args=[job_id, pdb_hex, config],
                queue="refinement",
                link=_make_callback_signature(celery_app, original_audit_id, user_email)
            )

        else:
            raise ValueError(f"Unknown protocol: {protocol}")

        queue_status   = STATE_QUEUED
        celery_task_id = task.id

    except Exception as e:
        logger.warning(f"[{job_id}] Celery dispatch failed: {e}")
        queue_status   = "beta_pending"
        celery_task_id = None

    return {
        "job_id":            job_id,
        "protocol":          protocol,
        "rationale":         selection["rationale"],
        "estimated_minutes": estimated_minutes,
        "state":             queue_status,
        "celery_task_id":    celery_task_id,
        "original_audit_id": original_audit_id,
        "message": (
            f"Job queued. Poll /refinement/status/{job_id} for updates."
            if queue_status == STATE_QUEUED
            else "B2 GPU worker not yet deployed. Use B1 callback flow."
        )
    }


def _make_callback_signature(celery_app, original_audit_id, user_email):
    """Build Celery chord callback for auto re-audit after execution."""
    return celery_app.signature(
        "worker.tasks.auto_callback",
        kwargs={
            "original_audit_id": original_audit_id,
            "user_email":        user_email
        }
    )


def _generate_rosetta_xml(audit_id: str, failing_laws: list, weights: Dict) -> str:
    """Generate Rosetta FastRelax XML from scoreterm weights."""
    failing_str = ", ".join(sorted(failing_laws))
    return f"""<ROSETTASCRIPTS>
  <!--
    Toscanini FastRelax Protocol (Phase B2 — Auto-Generated)
    Audit ID  : {audit_id}
    Violations: {failing_str}
    Generated by: execution_engine.py
  -->
  <SCOREFXNS>
    <ScoreFunction name="toscanini_relax" weights="ref2015">
      <Reweight scoretype="rama_prepro"  weight="{weights['rama_weight']}"/>
      <Reweight scoretype="fa_dun"       weight="{weights['dun_weight']}"/>
      <Reweight scoretype="fa_rep"       weight="{weights['rep_weight']}"/>
      <Reweight scoretype="omega"        weight="{weights['omega_weight']}"/>
      <Reweight scoretype="dslf_fa13"    weight="{weights['dslf_weight']}"/>
      <Reweight scoretype="cart_bonded"  weight="{weights['geom_weight']}"/>
    </ScoreFunction>
  </SCOREFXNS>
  <MOVERS>
    <FastRelax name="relax" scorefxn="toscanini_relax" repeats="5"
               cartesian="false" delete_virtual_residues_after_FastRelax="true"/>
    <MinMover  name="minimize" scorefxn="toscanini_relax"
               type="lbfgs_armijo_nonmonotone" tolerance="0.001" max_iter="200"/>
  </MOVERS>
  <PROTOCOLS>
    <Add mover="relax"/>
    <Add mover="minimize"/>
  </PROTOCOLS>
  <OUTPUT scorefxn="toscanini_relax"/>
</ROSETTASCRIPTS>"""
