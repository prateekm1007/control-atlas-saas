#!/usr/bin/env python3
"""Fetch the LIVE pull-request body (fail closed; Round-25 auditor direction 2).

The Round-24 workflow tried to read the live PR body but fell back to the
EVENT PAYLOAD when the API call failed or returned an empty body. That is a
fail-OPEN path: a stale event payload is not proof of the current live
declaration, and the exact-head binding does not prove that the body was
retrieved successfully at all.

This script is the whole retrieval path. It:
  * reads the live body from the authenticated GitHub API;
  * FAILS CLOSED when the API cannot be read, when the response is malformed,
    when the body is missing, empty, or not a string;
  * writes the retrieved body to --out and prints only fixed metadata
    (the byte count) -- never the body, never a header, never a credential;
  * supports --api-json for deterministic offline testing (the same pattern as
    the other checkers in this package).

Exit: 0 = the live body was retrieved and written; 1 = fail closed;
2 = usage error.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.error
import urllib.request

API = "https://api.github.com"


def fetch_pr(repo: str, number: str, token: str) -> dict:
    """The pull request record. Raises RuntimeError with a FIXED reason on any
    failure (no URL, no header, no credential, no response body)."""
    url = f"{API}/repos/{repo}/pulls/{number}"
    req = urllib.request.Request(
        url, method="GET",
        headers={"Accept": "application/vnd.github+json",
                 "X-GitHub-Api-Version": "2022-11-28",
                 **({"Authorization": "Bearer " + token} if token else {})})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw = resp.read()
    except urllib.error.HTTPError as e:
        # the status code only: no header, no body, no URL
        raise RuntimeError(f"the live PR API returned HTTP {e.code}") from None
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise RuntimeError(f"the live PR API could not be read ({type(e).__name__})") from None
    try:
        data = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        raise RuntimeError("the live PR API response is not parseable JSON") from None
    if not isinstance(data, dict):
        raise RuntimeError("the live PR API response is not a pull-request object")
    return data


def extract_body(pr) -> str:
    """The live body. A response that is not a pull-request object, an absent,
    non-string or empty body is a FAILURE."""
    if not isinstance(pr, dict):
        raise RuntimeError("the live PR API response is not a pull-request object")
    if "body" not in pr:
        raise RuntimeError("the live PR API response carries no body field")
    body = pr.get("body")
    if body is None:
        raise RuntimeError("the live pull-request body is null")
    if not isinstance(body, str):
        raise RuntimeError("the live pull-request body is not a string")
    if not body.strip():
        raise RuntimeError("the live pull-request body is empty")
    return body


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", required=True, help="owner/name of the repository")
    ap.add_argument("--pr", required=True, help="the pull request number")
    ap.add_argument("--out", required=True, help="file to write the live body to")
    ap.add_argument("--api-token-env", default="GH_TOKEN", metavar="VAR",
                    help="environment variable carrying the API token (default "
                         "GH_TOKEN; the value is never printed or written anywhere)")
    ap.add_argument("--api-json", help="offline fixture: a canned PR API response "
                                        "(tests); 'null' models an unreadable API")
    args = ap.parse_args()
    if not args.repo.strip() or not args.pr.strip():
        print("FAIL: --repo and --pr are required (fail closed)")
        return 2
    try:
        if args.api_json is not None:
            if args.api_json.lower() in ("null", "none", "error"):
                raise RuntimeError("the live PR API could not be read (fixture)")
            with open(args.api_json, encoding="utf-8") as fh:
                pr = json.load(fh)
        else:
            token = os.environ.get(args.api_token_env, "") if args.api_token_env else ""
            if not token:
                print(f"FAIL: no API token in ${args.api_token_env}; the live PR body "
                      f"cannot be read and the claim gate fails closed")
                return 1
            pr = fetch_pr(args.repo, args.pr.strip(), token)
        body = extract_body(pr)
    except RuntimeError as e:
        print(f"FAIL: {e} -- the live pull-request body could not be retrieved, so the "
              f"claim gate fails closed (the event payload is NOT a substitute)")
        return 1
    except (OSError, ValueError) as e:
        print(f"FAIL: the live pull-request body could not be retrieved ({type(e).__name__}) "
              f"-- the claim gate fails closed (the event payload is NOT a substitute)")
        return 1
    try:
        with open(args.out, "w", encoding="utf-8", newline="") as fh:
            fh.write(body)
    except OSError as e:
        print(f"FAIL: the live body could not be written to the output file "
              f"({type(e).__name__}) -- fail closed")
        return 1
    print(f"live PR body retrieved: {len(body.encode('utf-8'))} bytes written to "
          f"{os.path.basename(args.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
