"""
Runner provenance manifest — v2 (Round-12 direction 4, Round-13 delivery).

Covers:
  - build_provenance_manifest composes every required block, including
    the explicit calibration_executed=False record marker (a provenance
    record must never be mistakable for a calibration run).
  - Package-install audit: a repository that merely CONTAINS an
    installed/untracked tos package is NOT the engine source and is
    never reported (work-tree containment + git-tracking ownership).
  - Installed-package route: distribution metadata recorded, code_sha
    None; a commit is never inferred from a containing checkout.
  - No Git metadata: recorded explicitly as unavailable; strict
    (authoritative) mode fails closed with SystemExit.
  - Clean/dirty tree states are recorded, never refused, by the record
    path; the strict path is exercised separately.
  - Git remote URLs are sanitized of credentials in HTTPS and scp-style
    forms (EC section 3); no secret values enter manifests or exception
    text.
"""
import json
import subprocess
from pathlib import Path
from urllib.parse import urlsplit

import pytest

import tos as tos_pkg
from calibration.runner import (
    _engine_repo_root,
    _engine_provenance,
    _git_context,
    _sanitize_git_url,
    build_provenance_manifest,
    main,
    write_provenance_manifest,
)


def _run_git(cwd, *args):
    return subprocess.run(["git", "-C", str(cwd), *args],
                          capture_output=True, text=True, check=True)


def _make_repo(path):
    path.mkdir(parents=True, exist_ok=True)
    _run_git(path, "init", "-q")
    _run_git(path, "config", "user.email", "t@t")
    _run_git(path, "config", "user.name", "t")
    (path / "README.md").write_text("x")
    _run_git(path, "add", "README.md")
    _run_git(path, "commit", "-qm", "init")
    return path


# ── credential sanitization (HTTPS + scp-style) ──────────────────────────

def test_sanitize_strips_userinfo():
    assert (_sanitize_git_url("https://user:secret@github.com/o/r.git")
            == "https://github.com/o/r.git")
    assert (_sanitize_git_url("https://token@github.com/o/r.git")
            == "https://github.com/o/r.git")
    assert (_sanitize_git_url("https://github.com/o/r.git")
            == "https://github.com/o/r.git")
    assert _sanitize_git_url("git@github.com:o/r.git") == "github.com:o/r.git"
    assert _sanitize_git_url("ssh://git@github.com/o/r.git") == \
        "ssh://github.com/o/r.git"


def test_git_context_sanitizes_poisoned_remotes(tmp_path, monkeypatch):
    repo = _make_repo(tmp_path / "repo")
    _run_git(repo, "remote", "add", "origin",
             "https://user:secret@github.com/o/r.git")
    _run_git(repo, "remote", "add", "upstream", "git@github.com:o/up.git")
    import calibration.runner as runner
    monkeypatch.setattr(runner, "_engine_repo_root", lambda: repo)
    ctx = _git_context()
    assert ctx["remotes"]["origin"] == "https://github.com/o/r.git"
    assert ctx["remotes"]["upstream"] == "github.com:o/up.git"
    dumped = json.dumps(ctx)
    assert "user:secret" not in dumped and "git@" not in dumped


# ── normal checkout ──────────────────────────────────────────────────────

def test_engine_repo_root_in_normal_checkout():
    """The imported tos package is owned by this repository (tracked),
    so the git-checkout provenance route resolves."""
    root = _engine_repo_root()
    assert root is not None and (root / ".git").exists()


def test_manifest_has_required_blocks():
    m = build_provenance_manifest()
    for key in ("schema", "generated_utc", "runner_version",
                "calibration_executed", "record_type", "engine",
                "git_context", "runtime", "dependencies", "policy_binding"):
        assert key in m, f"missing block: {key}"
    assert m["schema"] == "toscanini-runner-provenance/2.0.0"
    assert m["calibration_executed"] is False
    assert m["record_type"] == "provenance-record-only"
    assert m["engine"]["source_provenance"] == "git-checkout"
    assert len(m["engine"]["code_sha"]) == 40
    assert m["engine"]["code_sha_source"]
    assert "python" in m["runtime"]
    assert "python_implementation" in m["runtime"]
    assert "platform" in m["runtime"]


def test_manifest_remotes_carry_no_userinfo():
    m = build_provenance_manifest()
    assert m["git_context"]["remotes"], "expected at least one remote"
    for name, url in m["git_context"]["remotes"].items():
        if "://" in url:
            parts = urlsplit(url)
            assert parts.username is None and parts.password is None, \
                f"userinfo survived sanitization in remote {name}: {url}"
        authority = url.split("://")[-1].split("/")[0]
        assert "@" not in authority, \
            f"userinfo survived sanitization in remote {name}: {url}"


def test_dependencies_recorded_explicitly():
    deps = build_provenance_manifest()["dependencies"]
    assert "numpy" in deps and "gemmi" in deps
    for name, version in deps.items():
        assert version, f"dependency {name} recorded with empty version"


# ── dirty tree is recorded, not refused (record path) ────────────────────

def test_dirty_tree_recorded_not_refused(tmp_path, monkeypatch):
    repo = _make_repo(tmp_path / "repo")
    (repo / "uncommitted.txt").write_text("dirty")
    import calibration.runner as runner
    monkeypatch.setattr(runner, "_engine_repo_root", lambda: repo)
    prov = _engine_provenance(strict=False)
    assert prov["source_provenance"] == "git-checkout"
    assert prov["working_tree_clean"] is False
    assert any("uncommitted.txt" in line for line in prov["tree_status"])
    # the record CLI path still succeeds on a dirty tree
    out = tmp_path / "prov.json"
    assert main(["--emit-provenance", str(out)]) == 0


# ── package-install audit: containing repo is NOT the owner ──────────────

def test_unrelated_parent_repo_is_rejected(tmp_path, monkeypatch):
    """A site-packages-style copy nested inside an UNRELATED git repo
    must not be attributed that repo's commit (Round-12 direction 4)."""
    repo = _make_repo(tmp_path / "unrelated-repo")
    fake_pkg = repo / "venv" / "lib" / "python3.12" / "site-packages" \
        / "tos" / "__init__.py"
    fake_pkg.parent.mkdir(parents=True)
    fake_pkg.write_text("")  # present on disk, NOT tracked by the repo
    monkeypatch.setattr(tos_pkg, "__file__", str(fake_pkg))
    assert _engine_repo_root() is None
    prov = _engine_provenance(strict=False)
    assert prov["code_sha"] is None
    assert prov["source_provenance"] in ("unavailable", "installed-package")
    # the unrelated repo's HEAD must NOT appear anywhere in the record
    unrelated_head = _run_git(repo, "rev-parse", "HEAD").stdout.strip()
    assert unrelated_head not in json.dumps(prov)


def test_no_git_metadata_records_unavailable(tmp_path, monkeypatch):
    fake_pkg = tmp_path / "site-packages" / "tos" / "__init__.py"
    fake_pkg.parent.mkdir(parents=True)
    fake_pkg.write_text("")
    monkeypatch.setattr(tos_pkg, "__file__", str(fake_pkg))
    prov = _engine_provenance(strict=False)
    assert prov["source_provenance"] in ("unavailable", "installed-package")
    assert prov["code_sha"] is None
    m = build_provenance_manifest()
    assert m["engine"]["code_sha"] is None
    assert m["git_context"]["available"] is False
    assert m["git_context"]["remotes"] == {}
    assert m["calibration_executed"] is False


def test_strict_engine_provenance_fails_closed(tmp_path, monkeypatch):
    fake_pkg = tmp_path / "site-packages" / "tos" / "__init__.py"
    fake_pkg.parent.mkdir(parents=True)
    fake_pkg.write_text("")
    monkeypatch.setattr(tos_pkg, "__file__", str(fake_pkg))
    with pytest.raises(SystemExit) as exc_info:
        _engine_provenance()  # strict default: authoritative measurement
    # the failure text carries no credential material
    assert "token" not in str(exc_info.value).lower()


def test_installed_package_route_records_distribution(tmp_path, monkeypatch):
    fake_pkg = tmp_path / "site-packages" / "tos" / "__init__.py"
    fake_pkg.parent.mkdir(parents=True)
    fake_pkg.write_text("")
    monkeypatch.setattr(tos_pkg, "__file__", str(fake_pkg))
    import calibration.runner as runner
    monkeypatch.setattr(
        runner, "_installed_distribution",
        lambda: {"name": "toscanini-engine", "version": "1.2.3"})
    prov = _engine_provenance(strict=False)
    assert prov["source_provenance"] == "installed-package"
    assert prov["code_sha"] is None
    assert prov["installed_distribution"]["version"] == "1.2.3"
    # strict mode still refuses: an installed artifact has no source commit
    with pytest.raises(SystemExit):
        _engine_provenance()


# ── record roundtrip and CLI ─────────────────────────────────────────────

def test_write_provenance_manifest_roundtrip(tmp_path):
    out = tmp_path / "nested" / "provenance.json"
    written = write_provenance_manifest(out)
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded == written
    assert loaded["schema"] == "toscanini-runner-provenance/2.0.0"
    assert loaded["calibration_executed"] is False
    assert loaded["record_type"] == "provenance-record-only"
    assert loaded["engine"]["code_sha"] == written["engine"]["code_sha"]


def test_emit_provenance_cli_writes_manifest(tmp_path):
    out = tmp_path / "prov" / "demo.json"
    rc = main(["--emit-provenance", str(out)])
    assert rc == 0
    m = json.loads(out.read_text(encoding="utf-8"))
    assert m["schema"].startswith("toscanini-runner-provenance/")
    # record-only nature is explicit: this is NOT a calibration run
    assert m["calibration_executed"] is False
    assert m["record_type"] == "provenance-record-only"
    assert len(m["engine"]["code_sha"]) == 40
    assert "working_tree_clean" in m["engine"]
    # no credentials anywhere in the written artifact
    text = out.read_text(encoding="utf-8")
    assert "ghp_" not in text and ":secret@" not in text and "git@" not in text
