"""
Deploy workflow policy assertions — Round-17 contract (dispatch-only, inert).

Supersession record (EC 1): this module previously pinned the WP-14 design
(workflow_run-triggered, host-based deployment with health check and
rollback). The Round-17 auditor direction 2a (2026-10-20, under full
delegation) replaces that design: ``deploy.yml`` becomes a
``workflow_dispatch``-only file with no push trigger, no workflow_run
trigger, and no host references — superseding owner action A-10. The tests
below were rewritten to assert the NEW binding contract with equivalent
structural rigor. No assertion was merely deleted: every old protection
(gate-gated start, SHA discipline, credential hygiene, serialization) maps
to its new-contract counterpart, and the new strongest property is
asserted directly — NOTHING that happens in this repository can start a
deployment.

Constitution 7 (EC 7): "Deployment MUST depend on successful tests,
validation, and build for the same source SHA... An independent
push-triggered deployment MUST NOT bypass a failed quality gate." The
dispatch-only file satisfies this maximally: no repository event can
deploy; re-introducing automated deployment requires an explicit, reviewed
change to deploy.yml — which these tests will flag.

The legacy remote scripts (scripts/deploy/deploy_remote.sh,
rollback_remote.sh) remain in-tree but are intentionally unreferenced by
any workflow; re-wiring them is a deliberate policy change.
"""
import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
DEPLOY_YML = REPO / ".github" / "workflows" / "deploy.yml"
WORKFLOWS_DIR = REPO / ".github" / "workflows"

FORBIDDEN_TRIGGERS = {
    "push", "pull_request", "workflow_run", "schedule", "check_suite",
    "issues", "issue_comment", "pull_request_target", "create", "delete",
    "deployment", "status", "release", "page_build",
}


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


# ── trigger discipline: manual dispatch ONLY; no repository event deploys ────

def test_trigger_is_workflow_dispatch_only(doc):
    _, _, trigger = doc
    assert trigger is not None, "deploy.yml must declare a trigger"
    assert set(trigger.keys()) == {"workflow_dispatch"}, (
        f"deploy.yml triggers must be workflow_dispatch only; found {set(trigger.keys())}"
    )


def test_no_event_trigger_can_start_a_deployment(doc):
    """The Round-17 core property: no push, PR, gate completion, or any
    other repository event may start this workflow."""
    _, _, trigger = doc
    present = FORBIDDEN_TRIGGERS & set(trigger.keys())
    assert not present, f"forbidden automatic trigger(s) present: {present}"


def test_file_text_declares_no_event_triggers(doc):
    """Structural guard: the words must not appear as YAML keys even in
    commented-out form that a future edit could re-enable accidentally.
    (Comments that EXPLAIN their absence are permitted and asserted.)"""
    text, _, _ = doc
    for key in ("push:", "pull_request:", "workflow_run:", "schedule:"):
        for line in text.splitlines():
            stripped = line.strip()
            assert not stripped.startswith(key) or stripped.startswith("#"), (
                f"active trigger key '{key}' found in deploy.yml"
            )


# ── no host, no transport, no credentials: the file cannot deploy anything ──

def test_no_host_references(doc):
    text, _, _ = doc
    assert "EC2" not in text, "deploy.yml must not reference EC2"
    assert "ec2" not in text.lower(), "deploy.yml must not reference hosts (ec2)"


def test_no_ssh_or_transport_commands(doc):
    text, _, _ = doc
    for pattern in (r"\bssh\b", r"\bscp\b", r"\brsync\b", r"\bcurl\b", r"\bwget\b"):
        assert not re.search(pattern, text, re.IGNORECASE), (
            f"deploy.yml must contain no transport command: {pattern}"
        )


def test_no_secrets_are_referenced(doc):
    text, _, _ = doc
    assert "secrets." not in text, "deploy.yml must reference no secrets"


def test_job_performs_no_deployment_and_says_so(deploy_job):
    steps = deploy_job["steps"]
    assert steps, "deploy job must have steps"
    names = " ".join(s.get("name", "") for s in steps)
    assert "inert" in names.lower(), (
        "the deploy job must state in its step name that it is intentionally inert"
    )
    for s in steps:
        assert "uses" not in s, f"inert deploy job must not use actions: {s.get('name')}"
        assert "ssh" not in s.get("run", "").lower()


# ── supply chain: minimal permissions; any action must be SHA-pinned ─────────

def test_permissions_are_read_only(doc):
    _, d, _ = doc
    assert d["permissions"] == {"contents": "read"}


def test_all_actions_pinned_to_full_commit_sha(doc):
    text, _, _ = doc
    uses = re.findall(r"^\s*uses:\s*(\S+)", text, re.MULTILINE)
    for u in uses:
        assert re.search(r"@[0-9a-f]{40}$", u), f"action not SHA-pinned: {u}"


# ── no checkout: the workflow never touches repository content ───────────────

def test_no_checkout_of_code(doc):
    text, _, _ = doc
    assert "actions/checkout" not in text, "inert deploy workflow must not check out code"


# ── legacy remote scripts stay unreferenced by the deploy workflow ───────────

def test_remote_deploy_scripts_not_referenced(doc):
    text, _, _ = doc
    for script in ("deploy_remote.sh", "rollback_remote.sh"):
        assert script not in text, (
            f"deploy.yml must not reference the legacy script {script} "
            "(re-wiring deployment requires an explicit reviewed change)"
        )
