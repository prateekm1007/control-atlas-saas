"""
WP-14 workflow-policy assertions: deployment must be gate-gated, SHA-exact,
health-verified, and rollback-capable.

Constitution 7.2/7.3: "Deployment MUST depend on successful tests, validation,
and build for the same source SHA... An independent push-triggered deployment
MUST NOT bypass a failed quality gate."

These tests parse .github/workflows/deploy.yml and the paired remote scripts
and assert the POLICY as structure — so any future edit that re-opens a
deploy bypass (failed gate, wrong branch, SHA drift, missing health check,
removed rollback) fails the CI gate itself.
"""
import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
DEPLOY_YML = REPO / ".github" / "workflows" / "deploy.yml"
DEPLOY_SH = REPO / "scripts" / "deploy" / "deploy_remote.sh"
ROLLBACK_SH = REPO / "scripts" / "deploy" / "rollback_remote.sh"


@pytest.fixture(scope="module")
def doc():
    text = DEPLOY_YML.read_text()
    d = yaml.safe_load(text)
    # PyYAML 1.1 parses the bare key `on:` as boolean True.
    trigger = d.get(True, d.get("on"))
    return text, d, trigger


@pytest.fixture(scope="module")
def deploy_job(doc):
    _, d, _ = doc
    return d["jobs"]["deploy"]


# ── trigger: only a completed main-branch gate run may start deployment ──────

def test_trigger_is_workflow_run_on_gate(doc):
    _, _, trigger = doc
    assert set(trigger.keys()) == {"workflow_run"}
    wr = trigger["workflow_run"]
    assert wr["workflows"] == ["Toscanini Stability Gate"]
    assert wr["types"] == ["completed"]
    assert wr["branches"] == ["main"]


def test_gate_failure_cannot_deploy(deploy_job):
    cond = deploy_job["if"]
    assert "github.event.workflow_run.conclusion == 'success'" in cond


def test_non_main_branch_cannot_deploy(deploy_job):
    cond = deploy_job["if"]
    assert "github.event.workflow_run.head_branch == 'main'" in cond


def test_workflow_identity_is_validated(deploy_job):
    cond = deploy_job["if"]
    assert "github.event.workflow_run.name == 'Toscanini Stability Gate'" in cond
    assert "github.event.workflow_run.path == '.github/workflows/toscanini-gate.yml'" in cond


# ── supply chain: minimal permissions, pinned action, no credentials to PRs ──

def test_permissions_are_read_only(doc):
    _, d, _ = doc
    assert d["permissions"] == {"contents": "read"}


def test_all_actions_pinned_to_full_commit_sha(doc):
    text, _, _ = doc
    uses = re.findall(r"uses:\s*(\S+)", text)
    assert uses, "no actions found"
    for u in uses:
        assert re.search(r"@[0-9a-f]{40}$", u), f"action not SHA-pinned: {u}"


def test_checkout_uses_gated_sha_and_no_persisted_credentials(deploy_job):
    steps = deploy_job["steps"]
    co = [s for s in steps if s.get("uses", "").startswith("actions/checkout")]
    assert len(co) == 1
    with_ref = co[0].get("with", {})
    assert with_ref.get("ref") == "${{ github.event.workflow_run.head_sha }}"
    assert with_ref.get("persist-credentials") is False


# ── SHA discipline: deploy the gated SHA, never drifting origin/main ────────

def test_no_reset_to_origin_main_anywhere(doc):
    text, _, _ = doc
    assert "reset --hard origin/main" not in text


def test_integrity_step_asserts_checkout_equals_gated_sha(deploy_job):
    steps = deploy_job["steps"]
    integ = [s for s in steps if "Integrity check" in s.get("name", "")]
    assert len(integ) == 1
    assert "head_sha\" ] != \"${DEPLOY_SHA}\"" in integ[0]["run"] or \
           "!= \"${DEPLOY_SHA}\"" in integ[0]["run"]


def test_remote_script_pins_exact_sha():
    src = DEPLOY_SH.read_text()
    assert 'DEPLOY_SHA="${1:?FATAL: deploy SHA argument required}"' in src
    assert 'git checkout --force "${DEPLOY_SHA}"' in src
    # No CODE may reset to origin/main (comments mentioning the prohibition are fine).
    code_lines = [l for l in src.splitlines() if not l.lstrip().startswith("#")]
    assert not any("origin/main" in l for l in code_lines), "remote script drifts to origin/main"
    assert '!= "${DEPLOY_SHA}"' in src  # remote mismatch guard


# ── fail-closed secrets: missing deployment secrets abort before SSH ────────

def test_missing_deployment_secrets_abort(deploy_job):
    preflight = [s for s in deploy_job["steps"] if s.get("name", "").startswith("Preflight")]
    assert len(preflight) == 1
    run = preflight[0]["run"]
    for secret in ("EC2_HOST", "EC2_USER", "EC2_SSH_KEY"):
        assert secret in run
    assert "exit 1" in run


# ── health: bounded verification against /health with healthy-status assert ──

def test_bounded_health_verification(deploy_job):
    health = [s for s in deploy_job["steps"] if "health verification" in s.get("name", "")]
    assert len(health) == 1
    run = health[0]["run"]
    assert "/health" in run
    assert 'assert d.get("status")=="healthy"' in run
    assert "attempts=10" in run and "interval=15" in run  # bounded, not unbounded


# ── rollback: automated path exists and is wired to failure ─────────────────

def test_automated_rollback_on_failure(deploy_job):
    rb = [s for s in deploy_job["steps"] if s.get("run", "").find("rollback_remote.sh") != -1]
    assert len(rb) == 1
    assert rb[0].get("if") == "failure()"
    assert "Automated rollback" in rb[0].get("name", "")


def test_rollback_script_restores_recorded_sha():
    src = ROLLBACK_SH.read_text()
    assert ".pre_deploy_sha" in src
    assert 'git checkout --force "${PRE_DEPLOY_SHA}"' in src


# ── credentials hygiene: 0600 key file + unconditional cleanup ───────────────

def test_ssh_key_installed_0600_and_cleaned(deploy_job):
    steps = deploy_job["steps"]
    install = [s for s in steps if "Install SSH key" in s.get("name", "")]
    cleanup = [s for s in steps if s.get("name", "").startswith("Remove SSH key")]
    assert len(install) == 1 and "install -m 600" in install[0]["run"]
    assert len(cleanup) == 1 and cleanup[0].get("if") == "always()"
    assert "rm -f /tmp/deploy_key" in cleanup[0]["run"]


# ── serialization: one production deploy at a time ───────────────────────────

def test_concurrency_serializes_deployments(doc):
    _, d, _ = doc
    assert d["concurrency"]["group"] == "deploy-production"
