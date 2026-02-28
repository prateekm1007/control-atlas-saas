"""
Toscanini B.4 — Notebook Generator Tests
Run: python3 dashboard/test_notebook_generator.py
"""
import sys, os, json

passed = 0
failed = 0

def run(name, fn):
    global passed, failed
    try:
        fn()
        print(f"  ✓ {name}")
        passed += 1
    except Exception as e:
        print(f"  ✗ {name}: {e}")
        failed += 1

# ── Minimal audit fixture ──────────────────────────────────────────────────────
AUDIT = {
    "governance": {"audit_id": "TESTAB12", "timestamp_utc": "2024-01-01T00:00:00+00:00"},
    "verdict":    {"binary": "VETO", "deterministic_score": 40,
                   "coverage_pct": 82.0, "det_passed": 8, "det_total": 12},
    "tier1":      {"laws": [
        {"law_id": "LAW-125", "status": "FAIL", "method": "deterministic",
         "observed": 8.2, "threshold": 2.0, "operator": "<="},
        {"law_id": "LAW-130", "status": "FAIL", "method": "deterministic",
         "observed": 12,  "threshold": 0,   "operator": "=="},
        {"law_id": "LAW-182", "status": "PASS", "method": "deterministic",
         "observed": 0.8, "threshold": 0.6, "operator": ">="},
    ]},
    "provenance":      {"source": "clash_structure.pdb"},
    "characterization":{"total_residues": 142},
}
TOKEN = "test.jwt.token.abc123"

print("\n=== NOTEBOOK GENERATOR TESTS (B.4) ===\n")

def t_openmm_kaggle_valid_json():
    from notebook_generator import generate_notebook
    nb = generate_notebook(AUDIT, TOKEN, "kaggle", "openmm")
    assert nb["nbformat"] == 4
    assert len(nb["cells"]) > 3
    # Must be serializable
    json.dumps(nb)

def t_rosetta_colab_valid_json():
    from notebook_generator import generate_notebook
    nb = generate_notebook(AUDIT, TOKEN, "colab", "rosetta")
    assert nb["nbformat"] == 4
    json.dumps(nb)

def t_token_in_openmm_notebook():
    from notebook_generator import generate_notebook
    nb = generate_notebook(AUDIT, TOKEN, "kaggle", "openmm")
    all_source = " ".join(
        c["source"] for c in nb["cells"] if c["cell_type"] == "code"
    )
    assert TOKEN in all_source, "Callback token not embedded in notebook"

def t_token_in_rosetta_notebook():
    from notebook_generator import generate_notebook
    nb = generate_notebook(AUDIT, TOKEN, "colab", "rosetta")
    all_source = " ".join(
        c["source"] for c in nb["cells"] if c["cell_type"] == "code"
    )
    assert TOKEN in all_source

def t_audit_id_in_notebook():
    from notebook_generator import generate_notebook
    nb = generate_notebook(AUDIT, TOKEN, "kaggle", "openmm")
    nb_str = json.dumps(nb)
    assert "TESTAB12" in nb_str

def t_bytes_output():
    from notebook_generator import generate_notebook_bytes
    b = generate_notebook_bytes(AUDIT, TOKEN, "kaggle", "openmm")
    assert isinstance(b, bytes)
    assert len(b) > 100
    # Must be valid JSON
    json.loads(b)

def t_protocol_auto_select_rosetta():
    from notebook_generator import select_protocol_for_audit
    # LAW-125 and LAW-130 → rosetta
    assert select_protocol_for_audit(AUDIT) == "rosetta"

def t_protocol_auto_select_openmm():
    from notebook_generator import select_protocol_for_audit
    audit_openmm = {**AUDIT, "tier1": {"laws": [
        {"law_id": "LAW-182", "status": "FAIL", "method": "deterministic",
         "observed": 0.3, "threshold": 0.6, "operator": ">="},
    ]}}
    assert select_protocol_for_audit(audit_openmm) == "openmm"

def t_invalid_platform_raises():
    from notebook_generator import generate_notebook
    try:
        generate_notebook(AUDIT, TOKEN, "aws_sagemaker", "openmm")
        assert False, "Should have raised"
    except AssertionError:
        pass

def t_platform_info_structure():
    from notebook_generator import get_platform_info
    for p in ("kaggle", "colab", "paperspace"):
        info = get_platform_info(p)
        assert "name"       in info
        assert "gpu_type"   in info
        assert "free_hours" in info
        assert "signup_url" in info

def t_redirect_url_contains_platform():
    from notebook_generator import get_platform_redirect_url
    url = get_platform_redirect_url("colab", "test_nb.ipynb")
    assert "colab.research.google.com" in url
    assert "test_nb.ipynb" in url

def t_kaggle_redirect_url():
    from notebook_generator import get_platform_redirect_url
    url = get_platform_redirect_url("kaggle", "test_nb.ipynb")
    assert "kaggle.com" in url

def t_no_stochastic_in_openmm_cells():
    """Verify no random seeds in generated OpenMM code."""
    from notebook_generator import generate_notebook
    nb = generate_notebook(AUDIT, TOKEN, "kaggle", "openmm")
    all_code = " ".join(
        c["source"] for c in nb["cells"] if c["cell_type"] == "code"
    )
    forbidden = ["random.seed", "np.random.seed", "time.time()"]
    for f in forbidden:
        assert f not in all_code, f"Stochastic call found: {f}"

def t_toscanini_metadata_in_notebook():
    from notebook_generator import generate_notebook
    nb = generate_notebook(AUDIT, TOKEN, "paperspace", "openmm")
    assert "toscanini" in nb["metadata"]
    assert nb["metadata"]["toscanini"]["audit_id"] == "TESTAB12"
    assert nb["metadata"]["toscanini"]["protocol"] == "openmm"

run("OpenMM/Kaggle notebook is valid JSON",        t_openmm_kaggle_valid_json)
run("Rosetta/Colab notebook is valid JSON",        t_rosetta_colab_valid_json)
run("Callback token embedded in OpenMM notebook",  t_token_in_openmm_notebook)
run("Callback token embedded in Rosetta notebook", t_token_in_rosetta_notebook)
run("Audit ID present in notebook",                t_audit_id_in_notebook)
run("generate_notebook_bytes returns bytes",        t_bytes_output)
run("Auto-select rosetta for clash/rama violations",t_protocol_auto_select_rosetta)
run("Auto-select openmm for burial violations",    t_protocol_auto_select_openmm)
run("Invalid platform raises AssertionError",      t_invalid_platform_raises)
run("Platform info has required fields",           t_platform_info_structure)
run("Colab redirect URL correct",                  t_redirect_url_contains_platform)
run("Kaggle redirect URL correct",                 t_kaggle_redirect_url)
run("No stochastic calls in OpenMM cells",         t_no_stochastic_in_openmm_cells)
run("Toscanini metadata in notebook",              t_toscanini_metadata_in_notebook)

print(f"\n=== RESULTS: {passed}/{passed+failed} passed ===\n")
if failed > 0:
    sys.exit(1)
