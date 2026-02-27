"""
dashboard/remediation_generator.py
Toscanini Remediation Script Generator (v2 - Deterministic)

Converts a Toscanini audit result into downloadable remediation artifacts.
No compute is performed. All outputs are scientifically parameterized
templates for execution in the user's own validated environment.

INVARIANTS:
- Deterministic output: same audit → byte-identical ZIP
- No execution of Rosetta, OpenMM, or any molecular dynamics
- No modification of law thresholds or governance logic
- All text is ASCII-safe
- Every artifact carries a mandatory recertification disclaimer

DETERMINISM GUARANTEES (v2):
- Fixed ZIP timestamps (2020-01-01 00:00:00)
- Sorted JSON keys in remediation_report.json
- Alphabetically sorted file insertion
- Conditional artifact inclusion (loop modeling only if LAW-110/160 fail)
- No random seeds, no stochastic flags
"""

import io
import json
import zipfile
from zipfile import ZipInfo
from datetime import datetime, timezone


# ── Deterministic ZIP metadata ───────────────────────────────────────────
FIXED_ZIP_DATE = (2020, 1, 1, 0, 0, 0)  # Normalized timestamp for all entries


# ── Law-to-remediation static mapping ────────────────────────────────────
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

Determinism: This artifact is generated deterministically. The same audit
input will always produce byte-identical output.
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


def _validate_no_stochastic_rosetta(xml_content):
    """Ensure Rosetta XML contains no random seeds or stochastic flags."""
    forbidden = ["random_seed", "constant_seed", "jran", "-run:constant_seed"]
    for flag in forbidden:
        if flag in xml_content:
            raise ValueError(f"STOCHASTIC FLAG DETECTED: '{flag}' in Rosetta XML. Artifact generation aborted.")


def _validate_openmm_deterministic(script_content):
    """Ensure OpenMM script is deterministic."""
    forbidden = ["random.seed", "np.random", "time.time()", "datetime.now"]
    for pattern in forbidden:
        if pattern in script_content:
            raise ValueError(f"NON-DETERMINISTIC CALL DETECTED: '{pattern}' in OpenMM script. Artifact generation aborted.")


def generate_rosetta_xml(audit_result):
    """Generate FastRelax XML. Returns string."""
    meta       = _get_audit_meta(audit_result)
    failing    = _get_failing_laws(audit_result)
    fail_ids   = [l["law_id"] for l in failing]

    rama_weight   = "2.0" if "LAW-125" in fail_ids else "1.0"
    dun_weight    = "2.0" if "LAW-150" in fail_ids else "1.0"
    rep_weight    = "1.5" if "LAW-130" in fail_ids else "0.55"
    omega_weight  = "2.0" if "LAW-135" in fail_ids else "1.0"
    geom_weight   = "2.0" if any(l in fail_ids for l in ["LAW-100", "LAW-120"]) else "1.0"
    dslf_weight   = "2.0" if "LAW-195" in fail_ids else "1.0"
    relax_cycles  = "5"
    failing_str   = ", ".join(fail_ids) if fail_ids else "NONE"

    # Extract residue-level diagnostics from LAW-125 (even if PASS)
    rama_residues = []
    all_laws = audit_result.get("tier1", {}).get("laws", [])
    for law in all_laws:
        if law.get("law_id") == "LAW-125" and law.get("granularity") == "residue":
            rama_residues = law.get("failing_residues", [])
            break
    
    rama_comment = ""
    if rama_residues:
        display_residues = rama_residues[:20]  # max 20 for XML clarity
        rama_str = ", ".join(map(str, display_residues))
        if len(rama_residues) > 20:
            rama_str += f" (+ {len(rama_residues) - 20} more)"
        rama_comment = f"    Ramachandran outliers at residues: {rama_str}\n"

    # Extract residue-level diagnostics from LAW-150
    rotamer_residues = []
    for law in all_laws:
        if law.get("law_id") == "LAW-150" and law.get("granularity") == "residue":
            rotamer_residues = law.get("failing_residues", [])
            break

    rotamer_comment = ""
    if rotamer_residues:
        rot_display = rotamer_residues[:20]
        rotamer_str = ", ".join(map(str, rot_display))
        if len(rotamer_residues) > 20:
            rotamer_str += f" (+ {len(rotamer_residues) - 20} more)"
        rotamer_comment = f"    Rotamer outliers at residues: {rotamer_str}\n"

    # Extract residue-level diagnostics from LAW-130 (clash pairs)
    clash_pairs = []
    for law in all_laws:
        if law.get("law_id") == "LAW-130" and law.get("granularity") == "residue_pair":
            clash_pairs = law.get("failing_residue_pairs", [])
            break

    clash_comment = ""
    if clash_pairs:
        display_pairs = clash_pairs[:10]  # max 10 pairs for XML clarity
        clash_str = ", ".join([f"{p[0]}-{p[1]}" for p in display_pairs])
        if len(clash_pairs) > 10:
            clash_str += f" (+ {len(clash_pairs) - 10} more pairs)"
        clash_comment = f"    Steric clashes between residue pairs: {clash_str}\n"

    # Extract residue-level diagnostics from LAW-135 (Omega planarity)
    omega_residues = []
    for law in all_laws:
        if law.get("law_id") == "LAW-135" and law.get("granularity") == "residue":
            omega_residues = law.get("failing_residues", [])
            break

    omega_comment = ""
    if omega_residues:
        display_residues = omega_residues[:15]  # max 15 for XML clarity
        omega_str = ", ".join(map(str, display_residues))
        if len(omega_residues) > 15:
            omega_str += f" (+ {len(omega_residues) - 15} more)"
        omega_comment = f"    Omega planarity outliers at residues: {omega_str}\n"
        omega_comment += f"    Recommendation: Constrain omega dihedrals during minimization\n"

    # Extract residue-level diagnostics from LAW-145 (Chirality)
    chiral_residues = []
    for law in all_laws:
        if law.get("law_id") == "LAW-145" and law.get("granularity") == "residue":
            chiral_residues = law.get("failing_residues", [])
            break

    chiral_comment = ""
    if chiral_residues:
        display_residues = chiral_residues[:10]  # max 10 for XML clarity
        chiral_str = ", ".join(map(str, display_residues))
        if len(chiral_residues) > 10:
            chiral_str += f" (+ {len(chiral_residues) - 10} more)"
        chiral_comment = f"    Chirality violations at residues: {chiral_str}\n"
        chiral_comment += f"    CRITICAL: Manual inspection required - automated refinement cannot fix chirality errors\n"

    xml = f"""<ROSETTASCRIPTS>
  <!--
    Toscanini FastRelax Protocol (Deterministic)
    Audit ID  : {meta['audit_id']}
    Generated : {meta['timestamp']}
    Verdict   : {meta['verdict']}
    Coverage  : {meta['coverage']}%
    Violations: {failing_str}
    {rama_comment}{rotamer_comment}{clash_comment}{omega_comment}{chiral_comment}    Residues  : {meta['n_residues']}
    Source    : {meta['source']}
{DISCLAIMER}
  -->

  <SCOREFXNS>
    <ScoreFunction name="toscanini_relax" weights="ref2015">
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

    _validate_no_stochastic_rosetta(xml)
    return xml


def generate_rosetta_flags(audit_result):
    """Generate Rosetta flags file. Returns string."""
    meta     = _get_audit_meta(audit_result)
    pdb_name = meta["source"].replace(" ", "_").replace("/", "_")

    flags = f"""# Toscanini FastRelax Flags (Deterministic)
# Audit ID  : {meta['audit_id']}
# Generated : {meta['timestamp']}
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
    """Generate OpenMM equilibration script. Returns string."""
    meta      = _get_audit_meta(audit_result)
    failing   = _get_failing_laws(audit_result)
    fail_ids  = [l["law_id"] for l in failing]
    failing_str = ", ".join(fail_ids) if fail_ids else "NONE"

    n_violations = len(_get_failing_det_laws(audit_result))
    sim_ns    = 2 if n_violations <= 2 else 5
    steps     = sim_ns * 500000
    pdb_name  = meta["source"].replace(" ", "_").replace("/", "_")

    script = f'''"""
Toscanini OpenMM Equilibration Script (Deterministic)
Audit ID  : {meta['audit_id']}
Generated : {meta['timestamp']}
Verdict   : {meta['verdict']}
Coverage  : {meta['coverage']}%
Violations: {failing_str}

{DISCLAIMER}

Requirements: OpenMM >= 7.7, PDBFixer >= 1.8
Usage: python openmm_equilibrate.py
"""

from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

# ── INPUT (DETERMINISTIC) ──────────────────────────────────────────────
PDB_FILE      = "{pdb_name}.pdb"
OUTPUT_PDB    = "{pdb_name}_equilibrated.pdb"
FORCEFIELD    = "amber14-all.xml"
WATER_MODEL   = "amber14/tip3pfb.xml"
TEMPERATURE_K = 300
TIMESTEP_FS   = 2
SIM_STEPS     = {steps}  # {sim_ns} ns at 2 fs/step
REPORT_EVERY  = 5000

print("Toscanini OpenMM Equilibration Protocol")
print("Audit ID: {meta['audit_id']}")
print("Violations targeted: {failing_str}")
print("Simulation: {{:.1f}} ns".format(SIM_STEPS * TIMESTEP_FS / 1e6))

# ── STRUCTURE PREPARATION ──────────────────────────────────────────────
pdb = PDBFile(PDB_FILE)
forcefield = ForceField(FORCEFIELD, WATER_MODEL)

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.0)
modeller.addSolvent(forcefield, model="tip3p", padding=10*angstroms)

# ── SYSTEM SETUP (DETERMINISTIC) ───────────────────────────────────────
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometers,
    constraints=HBonds,
)

integrator = LangevinMiddleIntegrator(
    TEMPERATURE_K * kelvin,
    1.0 / picosecond,
    TIMESTEP_FS * femtoseconds,
)

# ── SIMULATION ─────────────────────────────────────────────────────────
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

# ── OUTPUT ─────────────────────────────────────────────────────────────
positions = simulation.context.getState(getPositions=True).getPositions()
with open(OUTPUT_PDB, "w") as f:
    PDBFile.writeFile(simulation.topology, positions, f)

print(f"Refined structure saved to: {{OUTPUT_PDB}}")
print("NEXT STEP: Re-upload to Toscanini for independent re-certification.")
'''

    _validate_openmm_deterministic(script)
    return script


def generate_loop_modeling_script(audit_result):
    """Generate Rosetta KIC loop modeling XML. Returns string."""
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
    <Index name="loop_region" resnums="LOOP_START-LOOP_END"/>
  </RESIDUE_SELECTORS>

  <MOVERS>
    <LoopModeler name="kic_modeler"
                 scorefxn="score_loop"
                 protocol="kinematic"
                 loops_file="loops.txt">
    </LoopModeler>
  </MOVERS>

  <PROTOCOLS>
    <Add mover="kic_modeler"/>
  </PROTOCOLS>

</ROSETTASCRIPTS>"""
    return xml


def generate_remediation_report(audit_result):
    """Generate structured JSON remediation report. Returns dict."""
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

    fix_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "UNKNOWN": 3}
    prescriptions.sort(key=lambda x: (fix_order.get(x["fixability"], 3), x["law_id"]))

    report = {
        "toscanini_remediation_report": {
            "version":          "2.0",
            "generated_utc":    meta["timestamp"],
            "audit_id":         meta["audit_id"],
            "verdict":          meta["verdict"],
            "coverage_pct":     meta["coverage"],
            "deterministic_score": meta["det_score"],
            "total_residues":   meta["n_residues"],
            "source":           meta["source"],
            "total_violations": len(failing),
            "prescriptions":    prescriptions,
            "recertification_required": True,
            "determinism_guarantee": "Same audit input always produces byte-identical ZIP output",
            "disclaimer": (
                "This report is a prescription artifact. "
                "Toscanini performs no structural modification. "
                "All refined structures must be re-certified independently."
            ),
        }
    }
    return report


def generate_readme(audit_result, include_loop_modeling, callback_token=None):
    """Generate plain-text README. Returns string."""
    meta      = _get_audit_meta(audit_result)
    failing   = _get_failing_laws(audit_result)
    fail_ids  = [l["law_id"] for l in failing]
    pdb_name  = meta["source"].replace(" ", "_").replace("/", "_")

    loop_line = "  loop_modeling.xml        - KIC loop protocol (LAW-110/160 violations)"
    if not include_loop_modeling:
        loop_line = "  [loop_modeling.xml not included - no backbone gaps detected]"

    lines = [
        "TOSCANINI REMEDIATION PACKAGE (v2 - Deterministic)",
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
        "  README.txt               - This file",
        "  remediation_report.json  - Per-law prescriptions (sorted keys)",
        "  rosetta_relax.xml        - FastRelax protocol (parameterized)",
        "  rosetta.flags            - Rosetta execution flags",
        f"  {loop_line}",
        "  openmm_equilibrate.py    - OpenMM equilibration script",
        "",
        "DETERMINISM GUARANTEE",
        "-" * 40,
        "  This ZIP archive is generated deterministically.",
        "  The same audit input will always produce byte-identical output.",
        "  All timestamps are normalized to 2020-01-01 00:00:00.",
        "  JSON keys are sorted alphabetically.",
        "  Files are inserted in alphabetical order.",
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
        "       (only included if backbone gaps were detected)",
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

    if callback_token:
        sep = "=" * 60
        callback_section = (
            "\nMANAGED REFINEMENT CALLBACK (BETA)\n"
            + sep + "\n"
            + "  After executing refinement locally, re-upload your refined PDB\n"
            + "  using the callback token below.\n\n"
            + "  Callback Token (valid for 7 days):\n\n"
            + "    " + str(callback_token) + "\n\n"
            + "  Upload via Dashboard or curl.\n"
            + "  IMPORTANT: Token expires in 7 days. Single-use only.\n"
            + sep + "\n"
        )
        lines.insert(-5, callback_section)

    return "\n".join(lines)


def generate_remediation_zip(audit_result, callback_token=None):
    """
    Generate complete remediation package as deterministic ZIP.
    Same audit → byte-identical output.
    Returns raw bytes suitable for st.download_button().
    """
    meta    = _get_audit_meta(audit_result)
    failing = _get_failing_laws(audit_result)
    fail_ids = [l["law_id"] for l in failing]

    # Determine if loop modeling is needed
    needs_loop_modeling = any(lid in fail_ids for lid in ["LAW-110", "LAW-160"])

    # Build artifact list (will be sorted before insertion)
    artifacts = []

    # 1. Remediation report JSON (sorted keys for determinism)
    report = generate_remediation_report(audit_result)
    artifacts.append(("remediation_report.json", json.dumps(report, indent=2, sort_keys=True)))

    # 2. Rosetta FastRelax XML
    artifacts.append(("rosetta_relax.xml", generate_rosetta_xml(audit_result)))

    # 3. Rosetta flags
    artifacts.append(("rosetta.flags", generate_rosetta_flags(audit_result)))

    # 4. Loop modeling (conditional)
    if needs_loop_modeling:
        artifacts.append(("loop_modeling.xml", generate_loop_modeling_script(audit_result)))

    # 5. OpenMM script
    artifacts.append(("openmm_equilibrate.py", generate_openmm_script(audit_result)))

    # 6. README (must be generated last to know if loop modeling was included)
    artifacts.append(("README.txt", generate_readme(audit_result, needs_loop_modeling, callback_token)))

    # Sort alphabetically (deterministic across Python versions and dict insertion order)
    artifacts.sort(key=lambda x: x[0])

    # Write to ZIP with fixed timestamps
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for filename, content in artifacts:
            zi = ZipInfo(filename)
            zi.date_time = FIXED_ZIP_DATE
            zi.compress_type = zipfile.ZIP_DEFLATED
            zf.writestr(zi, content)

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
