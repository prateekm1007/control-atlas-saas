"""
Batch Processor: Multi-Structure Adjudication
----------------------------------------------
Processes a ZIP of PDB files through the Toscanini engine.

This module handles:
  - ZIP extraction (in-memory, no disk writes)
  - Per-structure adjudication via _run_physics_sync
  - Result aggregation and summary statistics

It has ZERO knowledge of HTTP or FastAPI.
It takes bytes in and returns structured results out.
"""

from __future__ import annotations

import io
import logging
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional

logger = logging.getLogger("toscanini.batch")


@dataclass
class StructureResult:
    """Result for a single structure within a batch."""
    filename: str
    candidate_id: str
    success: bool
    payload: Optional[dict] = None
    error: Optional[str] = None


@dataclass
class BatchSummary:
    """Aggregate statistics across all structures in a batch."""
    total: int = 0
    passed: int = 0
    vetoed: int = 0
    indeterminate: int = 0
    errors: int = 0
    mean_deterministic_score: float = 0.0
    common_failing_laws: list = field(default_factory=list)


@dataclass
class BatchResult:
    """Complete batch output."""
    results: list[StructureResult] = field(default_factory=list)
    summary: BatchSummary = field(default_factory=BatchSummary)


def _extract_pdb_files(zip_bytes: bytes) -> list[tuple[str, bytes]]:
    """
    Extract PDB files from a ZIP archive.

    Returns list of (filename, content_bytes) tuples.
    Skips non-PDB files, directories, and __MACOSX artifacts.
    """
    try:
        zf = zipfile.ZipFile(io.BytesIO(zip_bytes))
    except zipfile.BadZipFile as e:
        raise ValueError(f"Invalid ZIP file: {e}") from e

    pdb_files = []
    for info in zf.infolist():
        # Skip directories
        if info.is_dir():
            continue

        name = info.filename

        # Skip macOS artifacts and hidden files
        if "__MACOSX" in name or name.startswith(".") or "/." in name:
            continue

        # Only process .pdb and .ent files
        lower = name.lower()
        if not (lower.endswith(".pdb") or lower.endswith(".ent")):
            continue

        content = zf.read(info.filename)
        if len(content) < 50:
            logger.warning(f"Skipping {name}: too small ({len(content)} bytes)")
            continue

        pdb_files.append((name, content))

    if not pdb_files:
        raise ValueError("ZIP contains no valid PDB files (.pdb or .ent)")

    return pdb_files


def _derive_candidate_id(filename: str) -> str:
    """Extract a candidate ID from a filename."""
    # "structures/4HHB.pdb" → "4HHB"
    import os
    base = os.path.basename(filename)
    name, _ = os.path.splitext(base)
    return name


def process_batch(
    zip_bytes: bytes,
    run_physics_fn,
    mode: str = "Audit",
    t3_category: str = "NONE",
    max_structures: int = 100,
) -> BatchResult:
    """
    Process a batch of PDB structures from a ZIP file.

    Args:
        zip_bytes:       Raw bytes of the ZIP archive.
        run_physics_fn:  The physics engine function (injected to avoid import cycles).
                         Signature: (content_bytes, candidate_id, mode, t3_category) -> dict
        mode:            Audit mode string.
        t3_category:     Architecture category.
        max_structures:  Safety cap on number of structures to process.

    Returns:
        BatchResult with per-structure results and aggregate summary.
    """
    pdb_files = _extract_pdb_files(zip_bytes)

    if len(pdb_files) > max_structures:
        raise ValueError(
            f"Batch contains {len(pdb_files)} structures, "
            f"exceeding limit of {max_structures}"
        )

    logger.info(f"Batch: processing {len(pdb_files)} structures")

    results = []
    verdicts = []
    det_scores = []
    all_failing_laws = []

    for filename, content in pdb_files:
        candidate_id = _derive_candidate_id(filename)

        try:
            payload = run_physics_fn(content, candidate_id, mode, t3_category)

            verdict = payload.get("verdict", {}).get("binary", "ERROR")
            det_score = payload.get("verdict", {}).get("deterministic_score", 0)

            results.append(StructureResult(
                filename=filename,
                candidate_id=candidate_id,
                success=True,
                payload=payload,
            ))

            verdicts.append(verdict)
            det_scores.append(det_score)

            # Collect failing deterministic laws
            for law in payload.get("tier1", {}).get("laws", []):
                if law.get("method") == "deterministic" and law.get("status") != "PASS":
                    all_failing_laws.append(law["law_id"])

        except Exception as e:
            logger.error(f"Batch: {filename} failed: {e}")
            results.append(StructureResult(
                filename=filename,
                candidate_id=candidate_id,
                success=False,
                error=str(e),
            ))

    # Build summary
    verdict_counts = Counter(verdicts)
    law_counts = Counter(all_failing_laws)

    summary = BatchSummary(
        total=len(results),
        passed=verdict_counts.get("PASS", 0),
        vetoed=verdict_counts.get("VETO", 0),
        indeterminate=verdict_counts.get("INDETERMINATE", 0),
        errors=sum(1 for r in results if not r.success),
        mean_deterministic_score=round(
            sum(det_scores) / max(len(det_scores), 1), 1
        ),
        common_failing_laws=[
            {"law_id": lid, "count": count}
            for lid, count in law_counts.most_common(5)
        ],
    )

    return BatchResult(results=results, summary=summary)
