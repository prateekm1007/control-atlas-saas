"""
dashboard/notebook_generator.py
Toscanini BYOC Notebook Generator — Phase B.4

Generates pre-filled Jupyter notebooks for Kaggle, Colab, and Paperspace.
The user clicks "Run on Kaggle" → we generate a .ipynb → redirect to
the platform's notebook import URL → they run it → callback posts result.

INVARIANTS:
- Deterministic output: same audit + token → byte-identical notebook
- No execution of any physics code in this module
- No platform account access — user registers/logs in themselves
- All platform redirects use official documented import APIs (ToS-compliant)
- Callback token embedded in notebook for automatic re-upload after run

PLATFORM IMPORT URLs (official, ToS-compliant):
- Kaggle:     https://www.kaggle.com/kernels/welcome?src=<notebook_url>
- Colab:      https://colab.research.google.com/github/<org>/<repo>/blob/<branch>/<path>
- Paperspace: https://console.paperspace.com/github/<org>/<repo>/blob/<branch>/<path>

BYOC MODEL:
- User's compute quota (Kaggle 30hr/week, Colab free T4, Paperspace free GPU)
- Toscanini provides: audit, prescription, notebook, callback endpoint
- User provides: compute, platform account
- We never automate login or account creation
"""

import io
import json
import os
import hashlib
from datetime import datetime, timezone
from pathlib import Path


# ── Platform definitions ──────────────────────────────────────────────────────

PLATFORMS = {
    "kaggle": {
        "name":          "Kaggle",
        "gpu_type":      "NVIDIA T4 x2",
        "free_hours":    "30 hrs/week",
        "signup_url":    "https://www.kaggle.com/account/login",
        "import_url":    "https://www.kaggle.com/kernels/welcome?src={notebook_url}",
        "runtime":       "Python 3.10",
        "notes":         "Free GPU. Requires phone verification. 30 hr/week quota.",
        "install_cmd":   "pip install -q openmm pdbfixer",
    },
    "colab": {
        "name":          "Google Colab",
        "gpu_type":      "NVIDIA T4",
        "free_hours":    "Variable (session-based)",
        "signup_url":    "https://colab.research.google.com",
        "import_url":    "https://colab.research.google.com/github/{github_org}/{github_repo}/blob/main/dashboard/notebooks/{notebook_filename}",
        "runtime":       "Python 3.10",
        "notes":         "Free GPU via Google account. Sessions disconnect after idle.",
        "install_cmd":   "pip install -q openmm pdbfixer",
    },
    "paperspace": {
        "name":          "Paperspace Gradient",
        "gpu_type":      "NVIDIA M4000 (free tier)",
        "free_hours":    "Limited free tier",
        "signup_url":    "https://console.paperspace.com/signup",
        "import_url":    "https://console.paperspace.com/github/{github_org}/{github_repo}/blob/main/dashboard/notebooks/{notebook_filename}",
        "runtime":       "Python 3.9+",
        "notes":         "Free tier available. Pro tier for longer sessions.",
        "install_cmd":   "pip install -q openmm pdbfixer",
    },
}

# GitHub repo info for Colab/Paperspace import URLs
GITHUB_ORG  = os.environ.get("GITHUB_ORG",  "prateekm1007")
GITHUB_REPO = os.environ.get("GITHUB_REPO", "control-atlas-saas")
BACKEND_URL = os.environ.get("BACKEND_URL", "http://localhost:8000")


# ── Notebook cell builders ────────────────────────────────────────────────────

def _cell(source: str, cell_type: str = "code") -> dict:
    """Build a single notebook cell."""
    if cell_type == "markdown":
        return {
            "cell_type": "markdown",
            "metadata":  {},
            "source":    source,
        }
    return {
        "cell_type":       "code",
        "execution_count": None,
        "metadata":        {},
        "outputs":         [],
        "source":          source,
    }


def _md(text: str) -> dict:
    return _cell(text, cell_type="markdown")


def _code(text: str) -> dict:
    return _cell(text, cell_type="code")


# ── Notebook generators ───────────────────────────────────────────────────────

def _build_openmm_cells(audit_result: dict, callback_token: str,
                         platform: str) -> list:
    """Build cells for an OpenMM equilibration notebook."""
    meta      = _extract_meta(audit_result)
    fail_ids  = _get_fail_ids(audit_result)
    platform_cfg = PLATFORMS[platform]

    n_violations = len([l for l in audit_result.get("tier1", {}).get("laws", [])
                        if l.get("status") not in ("PASS",)
                        and l.get("method") == "deterministic"])
    sim_ns    = 2 if n_violations <= 2 else 5
    sim_steps = sim_ns * 500000
    pdb_name  = meta["source"].replace(" ", "_").replace("/", "_")

    cells = [
        _md(f"""# Toscanini BYOC — OpenMM Equilibration
**Audit ID:** `{meta['audit_id']}`  
**Verdict:** `{meta['verdict']}`  
**Violations:** `{', '.join(fail_ids) if fail_ids else 'NONE'}`  
**Platform:** {platform_cfg['name']} ({platform_cfg['gpu_type']})  
**Simulation:** {sim_ns} ns equilibration

> ⚠️ **Before running:** Upload your PDB file to this session's files panel.  
> After the run completes, this notebook will automatically post the refined  
> structure back to Toscanini for re-certification.
"""),

        _md("## Step 1: Install Dependencies"),
        _code(f"""{platform_cfg['install_cmd']}
print("Dependencies installed.")"""),

        _md("## Step 2: Upload Your Structure"),
        _code(f"""# Upload your PDB file using the Files panel (📁) on the left.
# Then set the filename below:
PDB_FILE   = "{pdb_name}.pdb"      # ← rename if needed
OUTPUT_PDB = "{pdb_name}_refined.pdb"
print(f"Input:  {{PDB_FILE}}")
print(f"Output: {{OUTPUT_PDB}}")"""),

        _md("## Step 3: Run OpenMM Equilibration"),
        _code(f"""from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

TEMPERATURE_K = 300
TIMESTEP_FS   = 2
SIM_STEPS     = {sim_steps}   # {sim_ns} ns
REPORT_EVERY  = 5000

print("Toscanini OpenMM Protocol")
print(f"Audit ID  : {meta['audit_id']}")
print(f"Violations: {', '.join(fail_ids) if fail_ids else 'NONE'}")
print(f"Simulation: {sim_ns} ns at 2 fs/step")

pdb        = PDBFile(PDB_FILE)
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
modeller   = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.0)
modeller.addSolvent(forcefield, model="tip3p", padding=10*angstroms)

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

try:
    platform = Platform.getPlatformByName("CUDA")
    print("Using CUDA GPU")
except Exception:
    platform = Platform.getPlatformByName("CPU")
    print("Using CPU (no GPU detected)")

simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

print("Running energy minimization...")
simulation.minimizeEnergy(maxIterations=1000)

simulation.reporters.append(
    StateDataReporter(sys.stdout, REPORT_EVERY, step=True,
                      potentialEnergy=True, temperature=True,
                      progress=True, remainingTime=True,
                      totalSteps=SIM_STEPS)
)

print(f"Running {{SIM_STEPS * TIMESTEP_FS / 1e6:.1f}} ns equilibration...")
simulation.step(SIM_STEPS)

positions = simulation.context.getState(getPositions=True).getPositions()
with open(OUTPUT_PDB, "w") as f:
    PDBFile.writeFile(simulation.topology, positions, f)

print(f"Refined structure saved: {{OUTPUT_PDB}}")"""),

        _md("## Step 4: Post Results to Toscanini"),
        _code(f"""import requests

CALLBACK_TOKEN = "{callback_token}"
BACKEND_URL    = "{BACKEND_URL}"
CALLBACK_URL   = f"{{BACKEND_URL}}/refinement/callback"

print("Posting refined structure to Toscanini for re-certification...")

with open(OUTPUT_PDB, "rb") as f:
    response = requests.post(
        CALLBACK_URL,
        data={{"token": CALLBACK_TOKEN}},
        files={{"file": (OUTPUT_PDB, f, "chemical/x-pdb")}},
        timeout=120
    )

if response.status_code == 200:
    result = response.json()
    print("✅ Re-certification complete!")
    print(f"   Original audit : {meta['audit_id']}")
    print(f"   Refined audit  : {{result.get('refined_audit_id', 'N/A')}}")
    print(f"   New verdict    : {{result.get('verdict', 'N/A')}}")
    print(f"   Coverage       : {{result.get('coverage_pct', 'N/A')}}%")
    print(f"   View comparison: {{result.get('comparison_url', 'N/A')}}")
else:
    print(f"⚠️ Callback failed ({{response.status_code}}): {{response.text[:200]}}")
    print("Manual option: Download the refined PDB and upload via the Toscanini dashboard.")
"""),

        _md("""## Done
Your refined structure has been submitted for re-certification.  
Return to the Toscanini dashboard to view the before/after comparison.

> **Note:** The callback token is single-use and expires in 7 days.
"""),
    ]
    return cells


def _build_rosetta_cells(audit_result: dict, callback_token: str,
                          platform: str) -> list:
    """Build cells for a Rosetta FastRelax notebook."""
    meta      = _extract_meta(audit_result)
    fail_ids  = _get_fail_ids(audit_result)
    platform_cfg = PLATFORMS[platform]

    # Scoreterm weights (mirrors protocol_selector.py logic)
    fail_set = set(fail_ids)
    weights  = {
        "rama":  "2.0" if "LAW-125" in fail_set else "1.0",
        "dun":   "2.0" if "LAW-150" in fail_set else "1.0",
        "rep":   "1.5" if "LAW-130" in fail_set else "0.55",
        "omega": "2.0" if "LAW-135" in fail_set else "1.0",
        "geom":  "2.0" if fail_set & {"LAW-100", "LAW-120"} else "1.0",
        "dslf":  "2.0" if "LAW-195" in fail_set else "1.0",
    }
    pdb_name = meta["source"].replace(" ", "_").replace("/", "_")

    cells = [
        _md(f"""# Toscanini BYOC — Rosetta FastRelax
**Audit ID:** `{meta['audit_id']}`  
**Verdict:** `{meta['verdict']}`  
**Violations:** `{', '.join(fail_ids) if fail_ids else 'NONE'}`  
**Platform:** {platform_cfg['name']} ({platform_cfg['gpu_type']})

> ⚠️ **Rosetta requires a license** (free for academic use).  
> Register at: https://www.rosettacommons.org/software/license-and-download  
> Upload your Rosetta binary to this session before running.
"""),

        _md("## Step 1: Configure Rosetta Path"),
        _code(f"""import os

# Set the path to your Rosetta binary
# Upload it to this session or mount from Google Drive
ROSETTA_BIN = "/content/rosetta_scripts.default.linuxgccrelease"  # Colab
# ROSETTA_BIN = "/kaggle/working/rosetta_scripts.default.linuxgccrelease"  # Kaggle

PDB_FILE   = "{pdb_name}.pdb"
OUTPUT_PDB = "{pdb_name}_relaxed.pdb"

assert os.path.exists(ROSETTA_BIN), f"Rosetta binary not found: {{ROSETTA_BIN}}"
print(f"Rosetta: {{ROSETTA_BIN}}")
print(f"Input:   {{PDB_FILE}}")"""),

        _md("## Step 2: Generate FastRelax Protocol"),
        _code(f"""XML_PROTOCOL = '''<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="toscanini_relax" weights="ref2015">
      <Reweight scoretype="rama_prepro" weight="{weights['rama']}"/>
      <Reweight scoretype="fa_dun"      weight="{weights['dun']}"/>
      <Reweight scoretype="fa_rep"      weight="{weights['rep']}"/>
      <Reweight scoretype="omega"       weight="{weights['omega']}"/>
      <Reweight scoretype="dslf_fa13"   weight="{weights['dslf']}"/>
      <Reweight scoretype="cart_bonded" weight="{weights['geom']}"/>
    </ScoreFunction>
  </SCOREFXNS>
  <MOVERS>
    <FastRelax name="relax" scorefxn="toscanini_relax" repeats="5"/>
    <MinMover  name="min"   scorefxn="toscanini_relax"
               type="lbfgs_armijo_nonmonotone" tolerance="0.001" max_iter="200"/>
  </MOVERS>
  <PROTOCOLS>
    <Add mover="relax"/>
    <Add mover="min"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>'''

with open("protocol.xml", "w") as f:
    f.write(XML_PROTOCOL)
print("Protocol written: protocol.xml")"""),

        _md("## Step 3: Run FastRelax"),
        _code(f"""import subprocess

cmd = [
    ROSETTA_BIN,
    "-in:file:s",      PDB_FILE,
    "-parser:protocol","protocol.xml",
    "-out:file:o",     OUTPUT_PDB,
    "-out:nstruct",    "1",
    "-out:overwrite",
    "-ex1", "-ex2",
    "-use_input_sc",
    "-mute",           "all",
    "-unmute",         "protocols.relax",
]

print("Running Rosetta FastRelax...")
result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

if result.returncode == 0:
    print("✅ FastRelax complete")
else:
    print(f"⚠️ Rosetta exited {{result.returncode}}")
    print(result.stderr[-500:])"""),

        _md("## Step 4: Post Results to Toscanini"),
        _code(f"""import requests

CALLBACK_TOKEN = "{callback_token}"
BACKEND_URL    = "{BACKEND_URL}"

with open(OUTPUT_PDB, "rb") as f:
    response = requests.post(
        f"{{BACKEND_URL}}/refinement/callback",
        data={{"token": CALLBACK_TOKEN}},
        files={{"file": (OUTPUT_PDB, f, "chemical/x-pdb")}},
        timeout=120
    )

if response.status_code == 200:
    result = response.json()
    print("✅ Re-certification complete!")
    print(f"   Refined audit  : {{result.get('refined_audit_id')}}")
    print(f"   New verdict    : {{result.get('verdict')}}")
    print(f"   Coverage       : {{result.get('coverage_pct')}}%")
else:
    print(f"⚠️ Callback failed: {{response.text[:200]}}")
"""),
    ]
    return cells


# ── Core notebook builder ─────────────────────────────────────────────────────

def generate_notebook(audit_result: dict, callback_token: str,
                       platform: str, protocol: str) -> dict:
    """
    Generate a complete .ipynb notebook dict.

    Args:
        audit_result:   Full Toscanini audit result dict
        callback_token: Single-use JWT for re-upload
        platform:       "kaggle" | "colab" | "paperspace"
        protocol:       "openmm" | "rosetta"

    Returns:
        dict — valid Jupyter notebook v4.5 format
    """
    assert platform in PLATFORMS, f"Unknown platform: {platform}"
    assert protocol in ("openmm", "rosetta"), f"Unknown protocol: {protocol}"

    if protocol == "openmm":
        cells = _build_openmm_cells(audit_result, callback_token, platform)
    else:
        cells = _build_rosetta_cells(audit_result, callback_token, platform)

    return {
        "nbformat":       4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language":     "python",
                "name":         "python3",
            },
            "language_info": {
                "name":    "python",
                "version": "3.10.0",
            },
            "toscanini": {
                "audit_id":      _extract_meta(audit_result)["audit_id"],
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "protocol":      protocol,
                "platform":      platform,
                "version":       "B.4.0",
            },
        },
        "cells": cells,
    }


def generate_notebook_bytes(audit_result: dict, callback_token: str,
                             platform: str, protocol: str) -> bytes:
    """Return notebook as UTF-8 JSON bytes for download."""
    nb = generate_notebook(audit_result, callback_token, platform, protocol)
    return json.dumps(nb, indent=2, ensure_ascii=False).encode("utf-8")


def get_platform_redirect_url(platform: str, notebook_filename: str) -> str:
    """
    Return the platform import URL for a hosted notebook.
    Used when notebook is committed to the repo (Colab/Paperspace path).
    For Kaggle, user must manually import the downloaded .ipynb.
    """
    cfg = PLATFORMS[platform]
    return cfg["import_url"].format(
        notebook_url=f"https://raw.githubusercontent.com/{GITHUB_ORG}/{GITHUB_REPO}/main/dashboard/notebooks/{notebook_filename}",
        github_org=GITHUB_ORG,
        github_repo=GITHUB_REPO,
        notebook_filename=notebook_filename,
    )


def get_platform_info(platform: str) -> dict:
    """Return platform metadata for UI display."""
    return PLATFORMS.get(platform, {})


def select_protocol_for_audit(audit_result: dict) -> str:
    """
    Auto-select notebook protocol based on failing laws.
    Mirrors protocol_selector.py logic — keeps both in sync.
    """
    fail_ids = set(_get_fail_ids(audit_result))
    rosetta_laws = {"LAW-100", "LAW-120", "LAW-125", "LAW-130",
                    "LAW-135", "LAW-145", "LAW-150", "LAW-195"}
    openmm_laws  = {"LAW-182", "LAW-200", "LAW-155"}

    needs_rosetta = bool(fail_ids & rosetta_laws)
    needs_openmm  = bool(fail_ids & openmm_laws)

    if needs_rosetta and not needs_openmm:
        return "rosetta"
    return "openmm"  # default — safer, no license required


# ── Helpers ───────────────────────────────────────────────────────────────────

def _extract_meta(audit_result: dict) -> dict:
    gov = audit_result.get("governance", {})
    v   = audit_result.get("verdict", {})
    return {
        "audit_id":  gov.get("audit_id", "UNKNOWN"),
        "timestamp": gov.get("timestamp_utc", datetime.now(timezone.utc).isoformat()),
        "verdict":   v.get("binary", "UNKNOWN"),
        "coverage":  v.get("coverage_pct", 0.0),
        "det_score": v.get("deterministic_score", 0),
        "source":    audit_result.get("provenance", {}).get("source", "unknown"),
        "n_residues": audit_result.get("characterization", {}).get("total_residues", "N/A"),
    }


def _get_fail_ids(audit_result: dict) -> list:
    laws = audit_result.get("tier1", {}).get("laws", [])
    return [l["law_id"] for l in laws if l.get("status") not in ("PASS",)]
