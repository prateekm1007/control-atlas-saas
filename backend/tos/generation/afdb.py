"""
Metadata-first AlphaFold DB acquisition (P0 items 11, 12, 14).

Replaces the hardcoded `model_v4` file URLs, which returned HTTP 404 for
every target in the v24 run (the file bucket only serves the CURRENT
model version). Acquisition is now:

  1. GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot}
  2. select the entry whose `modelEntityId` equals the requested entity
     (default `AF-{uniprot}-F1`) — NEVER `data[0]`;
  3. among matches, take the entry with the highest `latestVersion`;
  4. download its `cifUrl`;
  5. detect the format from the CONTENT (mmCIF vs PDB); JSON is never
     accepted as a structure;
  6. return the bytes with the URL, entity, version, SHA-256, the
     field names present in the API response, and — when the requested
     entity no longer exists — the migrated identifier the API served
     instead (P0DTC2 now answers as AF-0000000365840314).

Field names were verified live against the API on 2026-10-19
(evidence/item-11/): modelEntityId, latestVersion, allVersions, cifUrl,
pdbUrl, bcifUrl, uniprotAccession, globalMetricValue are all present.
The direction's warning about old/new field coexistence is recorded:
`pdbUrl` still exists alongside `cifUrl`; this module codes against
`cifUrl` only, as directed.

Everything fails closed: any deviation raises AFDBUnavailable with a
structured record — no fabricated structures, no silent fallbacks.
"""
from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from typing import Callable, Optional, Tuple

import requests

logger = logging.getLogger("toscanini.afdb")

AFDB_API_TEMPLATE = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
UNIPROT_API_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{accession}.json"
DEFAULT_TIMEOUT_S = 15

# Item 12: the fields this module codes against, verified live.
REQUIRED_FIELDS = ("modelEntityId", "latestVersion", "cifUrl")


def _http_get(url: str, timeout: int = DEFAULT_TIMEOUT_S):
    """Default transport. Tests inject their own via ``http_get``."""
    return requests.get(url, timeout=timeout)


class AFDBUnavailable(Exception):
    """Fail-closed acquisition failure with a structured record."""

    def __init__(self, accession: str, reason: str, record: dict):
        self.accession = accession
        self.reason = reason
        self.record = record
        super().__init__(f"AlphaFold DB acquisition failed for {accession}: {reason}")


@dataclass(frozen=True)
class AcquiredPrediction:
    """A verified, digest-pinned AFDB structure with full provenance."""
    content: bytes
    fmt: str                      # "cif" or "pdb", from CONTENT
    model_entity_id: str
    latest_version: int
    source_url: str               # the exact URL the bytes came from
    sha256: str
    uniprot_accession: Optional[str] = None
    global_metric_value: Optional[float] = None
    fields_present: Tuple[str, ...] = ()
    requested_entity_id: str = ""


def detect_structure_format(content: bytes) -> str:
    """Detect PDB vs mmCIF from the CONTENT, never from a URL.

    JSON is explicitly rejected: a JSON API body is not a structure, and
    treating it as one is how fabricated structures used to reach the
    engine (the old dispatcher returned JSON bytes as a 'pdb' file
    whenever the API entry lacked ``pdbUrl``).
    """
    if not content or not content.strip():
        raise ValueError("empty content: not a structure")
    head = content.lstrip()[:512]
    if head.startswith(b"{") or head.startswith(b"["):
        raise ValueError(
            "JSON payload is not a structure format — refusing to score it")
    if head.startswith(b"data_") or head.startswith(b"#"):
        return "cif"
    if b"_atom_site." in content[:16384] or b"loop_" in head:
        return "cif"
    first_line = content.lstrip().splitlines()[0][:6]
    if first_line in (b"ATOM  ", b"HETATM", b"HEADER", b"TITLE ", b"REMARK",
                      b"COMPND", b"CRYST1", b"SEQRES", b"MODEL ", b"END   "):
        return "pdb"
    raise ValueError(
        f"unrecognized structure format (first bytes: {content[:24]!r})")


def normalize_accession(candidate_id: str) -> str:
    """Extract the bare UniProt accession from AF-xxx-F1 or bare forms."""
    clean = candidate_id.strip().upper()
    if clean.startswith("AF-"):
        clean = clean[3:]
    if "-F" in clean:
        head, _, tail = clean.rpartition("-F")
        if tail.isdigit():
            clean = head
    return clean


def fetch_prediction(candidate_id: str,
                     requested_entity_id: Optional[str] = None,
                     http_get: Optional[Callable] = None,
                     timeout: int = DEFAULT_TIMEOUT_S) -> AcquiredPrediction:
    """Metadata-first acquisition of one AlphaFold prediction.

    Raises AFDBUnavailable on every failure mode: API unreachable, empty
    response, missing entity (including the migrated-identifier case),
    missing cifUrl, download failure, or a payload that is not a
    structure (JSON refused).
    """
    get = http_get or _http_get
    accession = normalize_accession(candidate_id)
    entity = requested_entity_id or f"AF-{accession}-F1"
    record: dict = {"api_url": AFDB_API_TEMPLATE.format(accession=accession),
                    "requested_entity_id": entity}

    # 1. metadata lookup
    try:
        resp = get(record["api_url"], timeout)
    except Exception as exc:  # network/timeout
        raise AFDBUnavailable(accession, f"metadata request failed: {exc}", record)
    if getattr(resp, "status_code", 0) != 200:
        raise AFDBUnavailable(
            accession, f"metadata HTTP {getattr(resp, 'status_code', '?')}", record)
    try:
        entries = json.loads(resp.content.decode("utf-8", errors="strict"))
    except Exception as exc:
        raise AFDBUnavailable(accession, f"metadata is not valid JSON: {exc}", record)
    if not isinstance(entries, list) or not entries:
        raise AFDBUnavailable(accession, "no predictions returned", record)

    # 2. select by modelEntityId — never data[0]
    matches = [e for e in entries
               if isinstance(e, dict) and e.get("modelEntityId") == entity]
    if not matches:
        served = [e.get("modelEntityId") for e in entries if isinstance(e, dict)]
        record.update({
            "served_entity_ids": served,
            "migrated_identifier": served[0] if served else None,
            "note": ("requested entity no longer exists; the API serves "
                     "these identifiers instead (recorded, failing closed)"),
        })
        raise AFDBUnavailable(
            accession,
            f"requested entity {entity} not served (API returned: {served})",
            record)

    # 3. highest latestVersion among matches
    try:
        best = max(matches, key=lambda e: int(e.get("latestVersion", 0)))
    except (TypeError, ValueError):
        raise AFDBUnavailable(accession, "latestVersion not parseable", record)

    fields_present = tuple(sorted(k for k in best.keys()))
    record["selected_entity_id"] = best.get("modelEntityId")
    record["latest_version"] = best.get("latestVersion")
    record["fields_present"] = fields_present
    for missing in REQUIRED_FIELDS:
        if missing not in best:
            raise AFDBUnavailable(
                accession,
                f"API entry lacks required field '{missing}' "
                f"(fields present: {fields_present})", record)

    # 4. download cifUrl
    url = best["cifUrl"]
    record["cif_url"] = url
    try:
        dl = get(url, timeout)
    except Exception as exc:
        raise AFDBUnavailable(accession, f"structure download failed: {exc}", record)
    if getattr(dl, "status_code", 0) != 200:
        raise AFDBUnavailable(
            accession,
            f"structure download HTTP {getattr(dl, 'status_code', '?')}", record)
    content = dl.content
    if not content or len(content) < 200:
        raise AFDBUnavailable(accession, "structure download suspiciously small", record)

    # 5. format from content — JSON refused
    try:
        fmt = detect_structure_format(content)
    except ValueError as exc:
        raise AFDBUnavailable(accession, str(exc), record)

    acquired = AcquiredPrediction(
        content=content,
        fmt=fmt,
        model_entity_id=best["modelEntityId"],
        latest_version=int(best["latestVersion"]),
        source_url=url,
        sha256=hashlib.sha256(content).hexdigest(),
        uniprot_accession=best.get("uniprotAccession"),
        global_metric_value=best.get("globalMetricValue"),
        fields_present=fields_present,
        requested_entity_id=entity,
    )
    logger.info("AFDB acquisition: %s v%s via %s (%s, sha256 %s...)",
                acquired.model_entity_id, acquired.latest_version, url,
                fmt, acquired.sha256[:12])
    return acquired


def uniprot_entry_type(accession: str,
                       http_get: Optional[Callable] = None,
                       timeout: int = DEFAULT_TIMEOUT_S) -> Optional[str]:
    """Informational UniProt status check (item 14).

    Returns the entryType (e.g. 'UniProtKB reviewed' / 'Inactive') or
    None when the check itself fails. Informational only — acquisition
    gating is the AFDB API's job; this exists so cohort tooling can
    RECORD why an accession has no prediction (P62988: entryType
    'Inactive' — merged entry — verified live 2026-10-19).
    """
    get = http_get or _http_get
    url = UNIPROT_API_TEMPLATE.format(accession=accession.strip().upper())
    try:
        resp = get(url, timeout)
        if getattr(resp, "status_code", 0) != 200:
            return None
        data = json.loads(resp.content.decode("utf-8", errors="strict"))
        return data.get("entryType")
    except Exception:
        return None
