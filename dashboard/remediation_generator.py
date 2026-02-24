"""
dashboard/remediation_generator.py
Toscanini Remediation Script Generator

Converts a Toscanini audit result into downloadable remediation artifacts.
No compute is performed. All outputs are scientifically parameterized
templates for execution in the user's own validated environment.

INVARIANTS:
- No execution of Rosetta, OpenMM, or any molecular dynamics
- No modification of law thresholds or governance logic
- All outputs are deterministic given identical audit input
- All text is ASCII-safe
- Every artifact carries a mandatory recertification disclaimer
"""

import io
import json
import zipfile
from datetime import datetime, timezone


# ── Law-to-remediation static mapping ─────────────────────────────────────────
# Maps each law ID to: fixability tier, primary method, secondary method, and
# the specific Rosetta scoreterm most relevant to that violation type.
LAW_REMEDIATION_MAP = {
    "LAW-100": {
        "fixability": "HIGH",
        "primary":    "Rosetta FastRelax",
        "secondary":  "Geometry Idealization",
        "scoreterm":  "bond_geometry",
        "description": "Bond length outliers. FastRelax with geometry idealization corrects bond deviations.",
    },
    "LAW-110": {
        "fixability": "MEDIUM",
        "primary":    "Loop Modeling (KIC)",
        "secondary":  "Rosetta FastRelax",
        "scoreterm":  "chainbreak",
        "description": "Backbone gaps detected. True gaps require loop rebuilding; coordinate errors may resolve with relax.",
    },
    "LAW-120": {
        "fixability": "HIGH",
        "primary":    "Rosetta FastRelax",
        "secondary":  "Geometry Idealization",
        "scoreterm":  "angle_constraint",
        "description": "Bond angle RMSD exceeds threshold. FastRelax with idealization restores Engh-Huber geometry.",
    },
    "LAW-125": {
        "fixability": "HIGH",
        "primary":    "Rosetta FastRelax",
        "secondary":  "Short MD Equilibration",
        "scoreterm":  "rama_prepro",
        "description": "Ramachandran outliers detected. FastRelax with rama_prepro scoreterm drives phi/psi into allowed regions.",
    },
    "LAW-130": {
        "fixability": "HIGH",
        "primary":    "Rosetta FastRelax",
        "secondary":  "Short MD Equilibration",
        "scoreterm":  "fa_rep",
        "description": "Steric clashes detected. FastRelax minimizes fa_rep (Lennard-Jones repulsion) to resolve overlaps.",
    },
    "LAW-135": {
        "fixability": "HIGH",
        "primary":    "Rosetta FastRelax",
        "secondary":  "Geometry Idealization",
        "scoreterm":  "omega",
        "description": "Cis-peptide bonds outside allowed positions. FastRelax with omega constraint restores planarity.",
    },
    "LAW-145": {
        "fixability": "LOW",
        "primary":    "Manual Correction",
        "secondary":  "Coordinate Inspection",
        "scoreterm":  "N/A",
        "description": "D-amino acid chirality violations. These require manual inspection and cannot be resolved by automated relax.",
    },
    "LAW-150": {
        "fixability": "HIGH",
        "primary":    "Rosetta FastRelax",
        "secondary":  "Rotamer Optimization",
        "scoreterm":  "fa_dun",
        "description": "Rotamer outliers detected. FastRelax samples the Dunbrack rotamer library (fa_dun) to find low-energy conformers.",
    },
    "LAW-160": {
        "fixability": "MEDIUM",
        "primary":    "Loop Modeling (KIC)",
        "secondary":  "Rosetta FastRelax",
        "scoreterm":  "chainbreak",
        "description": "CA-CA distance violation. May indicate chain break or coordinate error. Inspect before automated fix.",
    },
    "LAW-170": {
        "fixability": "LOW",
        "primary":    "Manual Correction",
        "secondary":  "Residue Renaming",
        "scoreterm":  "N/A",
        "description": "Non-standard residues detected. Requires manual assessment of whether modification is intentional.",
    },
    "LAW-182": {
        "fixability": "MEDIUM",
        "primary":    "Short MD Equilibration",
        "secondary":  "Rosetta FastRelax",
        "scoreterm":  "fa_sol",
        "description": "Hydrophobic burial insufficient. MD equilibration allows core reorganization under thermal motion.",
    },
    "LAW-195": {
        "fixability": "HIGH",
        "primary":    "Rosetta FastRelax",
        "secondary":  "Disulfide Optimization",
        "scoreterm":  "dslf_fa13",
        "description": "Disulfide geometry deviation. FastRelax with dslf_fa13 scoreterm optimizes cysteine pair geometry.",
    },
    "LAW-155": {
        "fixability": "MEDIUM",
        "primary":    "Rosetta FastRelax",
        "secondary":  "Short MD Equilibration",
        "scoreterm":  "fa_atr",
        "description": "Voxel occupancy below threshold. Packing issues addressed by relax or MD equilibration.",
    },
    "LAW-200": {
        "fixability": "MEDIUM",
        "primary":    "Short MD Equilibration",
        "secondary":  "Rosetta FastRelax",
        "scoreterm":  "fa_atr",
        "description": "Packing quality exceeds threshold. MD or relax compacts the structure.",
    },
}

DISCLAIMER = """
================================================================================
TOSCANINI REMEDIATION ARTIFACT
================================================================================
Generated by: Toscanini Structural Governance Engine
Purpose: Computationally parameterized protocol template

MANDATORY DISCLAIMER:
This artifact is a PRESCRIPTION, not an execution. Toscanini does not perform
molecular dynamics, energy minimization, or any structural modification.

Before use:
  1. Validate this script in your institutional Rosetta/OpenMM environment
  2. Review all parameters against your structure's specific context
  3. After refinement, re-upload the refined structure to Toscanini
  4. A new independent certification is required for every refined structure
  5. Toscanini certification applies only to the exact coordinate file audited

Liability: The Toscanini engine makes no warranty regarding the suitability
of these parameters for any specific structural context.
================================================================================
"""


def _get_failing_laws(audit_result):
    """Extract all failing laws from audit result."""
    laws = audit_result.get("tier1", {}).get("laws", [])
    return [l for l in laws if l.get("status") not in ("PASS",)]


def _get_failing_det_laws(audit_result):
    """Extract only failing deterministic laws."""
    laws = audit_result.get("tier1", {}).get("laws", [])
    return [l for l in laws if l.get("status") not in ("PASS",) and l.get("method") == "deterministic"]


def _get_audit_meta(audit_result):
    """Extract core audit metadata for script headers."""
    gov = audit_result.get("governance", {})
    v   = audit_result.get("verdict", {})
    return {
        "audit_id":   gov.get("audit_id", "UNKNOWN"),
        "timestamp":  gov.get("timestamp_utc", datetime.now(timezone.utc).isoformat()),
        "verdict":    v.get("binary", "UNKNOWN"),
        "coverage":   v.get("coverage_pct", 0.0),
        "det_score":  v.get("deterministic_score", 0),
        "source":     audit_result.get("provenance", {}).get("source", "unknown"),
        "n_residues": audit_result.get("characterization", {}).get("total_residues", "N/A"),
    }


def generate_rosetta_xml(audit_result):
    """
    Generate a FastRelax XML protocol parameterized for this structure's violations.
    Scoreterm weights are elevated for violated law categories.
    Returns XML string.
    """
    meta       = _get_audit_meta(audit_result)
    failing    = _get_failing_laws(audit_result)
    fail_ids   = [l["law_id"] for l in failing]

    # Determine which scoreterms need emphasis based on violations
    rama_weight   = "2.0" if "LAW-125" in fail_ids else "1.0"
    dun_weight    = "2.0" if "LAW-150" in fail_ids else "1.0"
    rep_weight    = "1.5" if "LAW-130" in fail_ids else "0.55"
    omega_weight  = "2.0" if "LAW-135" in fail_ids else "1.0"
    geom_weight   = "2.0" if any(l in fail_ids for l in ["LAW-100", "LAW-120"]) else "1.0"
    dslf_weight   = "2.0" if "LAW-195" in fail_ids else "1.0"

    relax_cycles  = "5"
    failing_str   = ", ".join(fail_ids) if fail_ids else "NONE"

    xml = f"""<ROSETTASCRIPTS>
  <!--
    Toscanini FastRelax Protocol
    Audit ID  : {meta['audit_id']}
    Generated : {meta['timestamp']}
    Verdict   : {meta['verdict']}
    Coverage  : {meta['coverage']}%
    Violations: {failing_str}
    Residues  : {meta['n_residues']}
    Source    : {meta['source']}
{DISCLAIMER}
  -->

  <SCOREFXNS>
    <ScoreFunction name="toscanini_relax" weights="ref2015">
      <!--
        Scoreterm weights elevated for detected violations.
        rama_prepro : Ramachandran potential  (LAW-125 weight={rama_weight})
        fa_dun      : Dunbrack rotamer energy (LAW-150 weight={dun_weight})
        fa_rep      : Lennard-Jones repulsion (LAW-130 weight={rep_weight})
        omega       : Peptide planarity       (LAW-135 weight={omega_weight})
        dslf_fa13   : Disulfide geometry      (LAW-195 weight={dslf_weight})
        cart_bonded : Bond/angle geometry     (LAW-100/120 weight={geom_weight})
      -->
      <Reweight scoretype="rama_prepro"  weight="{rama_weight}"/>
      <Reweight scoretype="fa_dun"       weight="{dun_weight}"/>
      <Reweight scoretype="fa_rep"       weight="{rep_weight}"/>
      <Reweight scoretype="omega"        weight="{omega_weight}"/>
      <Reweight scoretype="dslf_fa13"   weight="{dslf_weight}"/>
      <Reweight scoretype="cart_bonded"  weight="{geom_weight}"/>
    </ScoreFunction>
  </SCOREFXNS>

  <MOVERS>
    <FastRelax name="toscanini_fast_relax"
               scorefxn="toscanini_relax"
               repeats="{relax_cycles}"
               cartesian="false"
               delete_virtual_residues_after_FastRelax="true">
      <!--
        repeats={relax_cycles}: Conservative default. Increase to 10-15 for severe violations.
        cartesian=false: Torsion-space relax. Set true if bond geometry violations are primary.
      -->
    </FastRelax>

    <MinMover name="final_minimize"
              scorefxn="toscanini_relax"
              type="lbfgs_armijo_nonmonotone"
              tolerance="0.001"
              max_iter="200">
    </MinMover>
  </MOVERS>

  <PROTOCOLS>
    <Add mover="toscanini_fast_relax"/>
    <Add mover="final_minimize"/>
  </PROTOCOLS>

  <OUTPUT scorefxn="toscanini_relax"/>

</ROSETTASCRIPTS>"""
    return xml


def generate_rosetta_flags(audit_result):
    """
    Generate a Rosetta flags file for the FastRelax run.
    Returns flags string.
    """
    meta     = _get_audit_meta(audit_result)
    pdb_name = meta["source"].replace(" ", "_").replace("/", "_")

    flags = f"""# Toscanini FastRelax Flags
# Audit ID  : {meta['audit_id']}
# Generated : {meta['timestamp']}
# Verdict   : {meta['verdict']}
# Coverage  : {meta['coverage']}%
#
# Usage: rosetta_scripts.default.linuxgccrelease @rosetta.flags
#
-in:file:s {pdb_name}.pdb
-parser:protocol rosetta_relax.xml
-out:file:o {pdb_name}_refined.pdb
-out:nstruct 1
-out:overwrite
-ex1
-ex2
-use_input_sc
-flip_HNQ
-no_optH false
-relax:constrain_relax_to_start_coords false
-relax:ramp_constraints false
-ignore_unrecognized_res false
-mute all
-unmute protocols.relax
"""
    return flags


def generate_openmm_script(audit_result):
    """
    Generate a runnable OpenMM equilibration script parameterized for this structure.
    Returns Python script string.
    """
    meta      = _get_audit_meta(audit_result)
    failing   = _get_failing_laws(audit_result)
    fail_ids  = [l["law_id"] for l in failing]
    failing_str = ", ".join(fail_ids) if fail_ids else "NONE"

    # Adjust simulation length based on violation severity
    n_violations = len(_get_failing_det_laws(audit_result))
    sim_ns    = 2 if n_violations <= 2 else 5
    steps     = sim_ns * 500000  # 2 fs timestep

    pdb_name  = meta["source"].replace(" ", "_").replace("/", "_")

    script = f'''"""
Toscanini OpenMM Equilibration Script
Audit ID  : {meta['audit_id']}
Generated : {meta['timestamp']}
Verdict   : {meta['verdict']}
Coverage  : {meta['coverage']}%
Violations: {failing_str}
Residues  : {meta['n_residues']}

{DISCLAIMER}

Requirements: OpenMM >= 7.7, PDBFixer >= 1.8
Usage: python openmm_equilibrate.py
"""

from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

# ── INPUT ──────────────────────────────────────────────────────────────────
PDB_FILE      = "{pdb_name}.pdb"
OUTPUT_PDB    = "{pdb_name}_equilibrated.pdb"
FORCEFIELD    = "amber14-all.xml"
WATER_MODEL   = "amber14/tip3pfb.xml"
TEMPERATURE_K = 300
TIMESTEP_FS   = 2
SIM_STEPS     = {steps}  # {sim_ns} ns at 2 fs/step
REPORT_EVERY  = 5000     # Log every 10 ps

print("Toscanini OpenMM Equilibration Protocol")
print("Audit ID: {meta['audit_id']}")
print("Violations targeted: {failing_str}")
print("Simulation: {{:.1f}} ns".format(SIM_STEPS * TIMESTEP_FS / 1e6))

# ── STRUCTURE PREPARATION ──────────────────────────────────────────────────
pdb = PDBFile(PDB_FILE)
forcefield = ForceField(FORCEFIELD, WATER_MODEL)

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.0)
modeller.addSolvent(forcefield, model="tip3p", padding=10*angstroms)

# ── SYSTEM SETUP ───────────────────────────────────────────────────────────
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometers,
    constraints=HBonds,
)

# Integrator: Langevin with conservative timestep
integrator = LangevinMiddleIntegrator(
    TEMPERATURE_K * kelvin,
    1.0 / picosecond,
    TIMESTEP_FS * femtoseconds,
)

# ── SIMULATION ─────────────────────────────────────────────────────────────
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

print("Running energy minimization...")
simulation.minimizeEnergy(maxIterations=1000)
print("Energy minimization complete.")

simulation.reporters.append(
    StateDataReporter(sys.stdout, REPORT_EVERY, step=True,
                      potentialEnergy=True, temperature=True, progress=True,
                      remainingTime=True, totalSteps=SIM_STEPS)
)

print(f"Running {{SIM_STEPS * TIMESTEP_FS / 1e6:.1f}} ns equilibration...")
simulation.step(SIM_STEPS)
print("Equilibration complete.")

# ── OUTPUT ─────────────────────────────────────────────────────────────────
positions = simulation.context.getState(getPositions=True).getPositions()
with open(OUTPUT_PDB, "w") as f:
    PDBFile.writeFile(simulation.topology, positions, f)

print(f"Refined structure saved to: {{OUTPUT_PDB}}")
print("NEXT STEP: Re-upload to Toscanini for independent re-certification.")
print("Certification is required for every refined coordinate file.")
'''
    return script


def generate_loop_modeling_script(audit_result):
    """
    Generate a Rosetta KIC loop modeling script for backbone gap violations.
    Only relevant when LAW-110 or LAW-160 fails.
    Returns XML string.
    """
    meta     = _get_audit_meta(audit_result)
    failing  = _get_failing_laws(audit_result)
    fail_ids = [l["law_id"] for l in failing]

    xml = f"""<ROSETTASCRIPTS>
  <!--
    Toscanini Loop Modeling Protocol (KIC)
    Audit ID  : {meta['audit_id']}
    Generated : {meta['timestamp']}
    Violations: {", ".join(fail_ids)}

    NOTE: KIC loop modeling requires manual identification of loop boundaries.
    Inspect your structure to identify which residues span the backbone gap
    (LAW-110/LAW-160 violations) before executing this protocol.

    Replace LOOP_START and LOOP_END with actual residue numbers.
{DISCLAIMER}
  -->

  <SCOREFXNS>
    <ScoreFunction name="score_loop" weights="ref2015"/>
  </SCOREFXNS>

  <RESIDUE_SELECTORS>
    <!--
      EDIT REQUIRED: Replace LOOP_START and LOOP_END with residue numbers
      from your structure that correspond to the backbone gap region.
    -->
    <Index name="loop_region" resnums="LOOP_START-LOOP_END"/>
  </RESIDUE_SELECTORS>

  <MOVERS>
    <LoopModeler name="kic_modeler"
                 scorefxn="score_loop"
                 protocol="kinematic"
                 loops_file="loops.txt">
      <!--
        loops.txt format:
        LOOP [start_residue] [end_residue] [cut_point] 0 1
        Example: LOOP 42 55 48 0 1
      -->
    </LoopModeler>
  </MOVERS>

  <PROTOCOLS>
    <Add mover="kic_modeler"/>
  </PROTOCOLS>

</ROSETTASCRIPTS>"""
    return xml


def generate_remediation_report(audit_result):
    """
    Generate a structured JSON remediation report with per-law prescriptions.
    Returns dict (caller serializes to JSON).
    """
    meta        = _get_audit_meta(audit_result)
    all_laws    = audit_result.get("tier1", {}).get("laws", [])
    failing     = _get_failing_laws(audit_result)

    prescriptions = []
    for law in failing:
        lid     = law.get("law_id", "UNKNOWN")
        remap   = LAW_REMEDIATION_MAP.get(lid, {
            "fixability":  "UNKNOWN",
            "primary":     "Manual Review",
            "secondary":   "N/A",
            "scoreterm":   "N/A",
            "description": "No automated prescription available for this law.",
        })
        prescriptions.append({
            "law_id":        lid,
            "title":         law.get("title", ""),
            "method":        law.get("method", ""),
            "observed":      law.get("observed", "N/A"),
            "threshold":     law.get("threshold", "N/A"),
            "operator":      law.get("operator", ""),
            "status":        law.get("status", ""),
            "fixability":    remap["fixability"],
            "primary_fix":   remap["primary"],
            "secondary_fix": remap["secondary"],
            "scoreterm":     remap["scoreterm"],
            "description":   remap["description"],
        })

    # Sort: HIGH fixability first, then MEDIUM, then LOW
    fix_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "UNKNOWN": 3}
    prescriptions.sort(key=lambda x: fix_order.get(x["fixability"], 3))

    report = {
        "toscanini_remediation_report": {
            "version":          "1.0",
            "generated_utc":    datetime.now(timezone.utc).isoformat(),
            "audit_id":         meta["audit_id"],
            "verdict":          meta["verdict"],
            "coverage_pct":     meta["coverage"],
            "deterministic_score": meta["det_score"],
            "total_residues":   meta["n_residues"],
            "source":           meta["source"],
            "total_violations": len(failing),
            "prescriptions":    prescriptions,
            "recertification_required": True,
            "disclaimer": (
                "This report is a prescription artifact. "
                "Toscanini performs no structural modification. "
                "All refined structures must be re-certified independently."
            ),
        }
    }
    return report


def generate_readme(audit_result):
    """Generate plain-text README for the remediation package."""
    meta      = _get_audit_meta(audit_result)
    failing   = _get_failing_laws(audit_result)
    fail_ids  = [l["law_id"] for l in failing]
    pdb_name  = meta["source"].replace(" ", "_").replace("/", "_")

    lines = [
        "TOSCANINI REMEDIATION PACKAGE",
        "=" * 60,
        "",
        f"Audit ID  : {meta['audit_id']}",
        f"Generated : {meta['timestamp']}",
        f"Verdict   : {meta['verdict']}",
        f"Coverage  : {meta['coverage']}%",
        f"Det. Score: {meta['det_score']}/100",
        f"Source    : {meta['source']}",
        f"Residues  : {meta['n_residues']}",
        f"Violations: {', '.join(fail_ids) if fail_ids else 'NONE'}",
        "",
        "PACKAGE CONTENTS",
        "-" * 40,
        "  remediation_report.json  - Per-law prescriptions and priorities",
        "  rosetta_relax.xml        - FastRelax protocol (parameterized)",
        "  rosetta.flags            - Rosetta execution flags",
        "  loop_modeling.xml        - KIC loop protocol (if LAW-110/160 fail)",
        "  openmm_equilibrate.py    - OpenMM equilibration script",
        "  README.txt               - This file",
        "",
        "EXECUTION ORDER",
        "-" * 40,
        "  1. Review remediation_report.json to understand all violations",
        "  2. For rotamer/clash/rama violations (HIGH fixability):",
        f"       rosetta_scripts.default.linuxgccrelease @rosetta.flags",
        "  3. For hydrophobic burial / packing violations (MEDIUM fixability):",
        f"       python openmm_equilibrate.py",
        "  4. For backbone gaps (LAW-110/LAW-160):",
        "       Edit loop_modeling.xml with correct residue numbers first",
        "  5. For chirality / non-standard residues (LOW fixability):",
        "       Manual inspection required. No automated fix available.",
        "",
        "MANDATORY RECERTIFICATION",
        "-" * 40,
        "  After any refinement step, re-upload the refined PDB to Toscanini.",
        "  Certification applies only to the exact coordinate file audited.",
        "  A new audit ID will be assigned to the refined structure.",
        "",
        "INPUT FILE NAMING",
        "-" * 40,
        f"  Rosetta expects: {pdb_name}.pdb",
        f"  OpenMM  expects: {pdb_name}.pdb",
        f"  Rename your PDB file to match, or edit the flags/script accordingly.",
        "",
        DISCLAIMER,
    ]
    return "\n".join(lines)


def generate_remediation_zip(audit_result):
    """
    Generate a complete remediation package as a ZIP archive.
    Returns raw bytes suitable for st.download_button().
    """
    meta    = _get_audit_meta(audit_result)
    failing = _get_failing_laws(audit_result)
    fail_ids = [l["law_id"] for l in failing]

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:

        # 1. Remediation report JSON
        report = generate_remediation_report(audit_result)
        zf.writestr(
            "remediation_report.json",
            json.dumps(report, indent=2),
        )

        # 2. Rosetta FastRelax XML
        zf.writestr("rosetta_relax.xml", generate_rosetta_xml(audit_result))

        # 3. Rosetta flags
        zf.writestr("rosetta.flags", generate_rosetta_flags(audit_result))

        # 4. Loop modeling (always include; user decides if relevant)
        zf.writestr("loop_modeling.xml", generate_loop_modeling_script(audit_result))

        # 5. OpenMM script
        zf.writestr("openmm_equilibrate.py", generate_openmm_script(audit_result))

        # 6. README
        zf.writestr("README.txt", generate_readme(audit_result))

    buf.seek(0)
    return buf.read()


def should_show_remediation(audit_result):
    """
    Returns True if this audit result warrants a remediation package.
    PASS + coverage >= 70 = certified, show dossier only.
    Anything else = show remediation package.
    """
    v   = audit_result.get("verdict", {})
    return not (
        v.get("binary") == "PASS" and v.get("coverage_pct", 0) >= 70
    )
