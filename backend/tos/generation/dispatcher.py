import requests
import logging
import re

from .afdb import fetch_prediction, AFDBUnavailable

logger = logging.getLogger("toscanini.dispatcher")


class SourceUnavailable(Exception):
    """WP-4: raised when every structure source for a candidate fails.

    The route layer turns this into HTTP 503 with the list of sources that
    were tried. No verdict is produced, nothing is written to the NKG, and
    no comparison is recorded — a fabricated structure must never be scored.
    """

    def __init__(self, candidate_id: str, sources: list):
        self.candidate_id = candidate_id
        self.sources = list(sources)
        super().__init__(
            f"All structure sources exhausted for '{candidate_id}' "
            f"(tried: {', '.join(self.sources)})"
        )


def _is_rcsb_id(candidate_id: str) -> bool:
    """Detect 4-character PDB IDs (e.g. 4HHB, 1CRN, 2PTC)."""
    clean = candidate_id.strip().upper()
    return bool(re.match(r'^[0-9][A-Z0-9]{3}$', clean))


def _is_af_id(candidate_id: str) -> bool:
    """Detect AlphaFold IDs (e.g. AF-P01308-F1 or bare UniProt P01308)."""
    clean = candidate_id.strip()
    if clean.startswith("AF-"):
        return True
    # Bare UniProt: 6-10 alphanumeric, starts with letter
    if re.match(r'^[A-Z][A-Z0-9]{4,9}$', clean.upper()) and not _is_rcsb_id(clean):
        return True
    return False


class GenerationDispatcher:
    """PILLAR 05: Industrial Acquisition — RCSB + AFDB dual-source downloader."""

    @staticmethod
    def acquire(candidate_id: str, fallback):
        clean = candidate_id.strip()

        if _is_rcsb_id(clean):
            return GenerationDispatcher._acquire_rcsb(clean.upper())
        else:
            return GenerationDispatcher._acquire_afdb(clean)

    @staticmethod
    def _acquire_rcsb(pdb_id: str):
        """Fetch experimental structure from RCSB PDB."""
        urls = [
            f"https://files.rcsb.org/download/{pdb_id}.pdb",
            f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb",
        ]

        for url in urls:
            try:
                logger.info(f"RCSB Acquisition: {url}")
                res = requests.get(url, timeout=30)
                if res.status_code == 200 and len(res.content) > 500:
                    atom_count = sum(1 for l in res.content.decode(errors="ignore").split("\n")
                                     if l.startswith("ATOM"))
                    logger.info(f"✅ RCSB Success: {pdb_id} — {atom_count} ATOM lines, {len(res.content)} bytes")
                    return (res.content, pdb_id, "pdb")
            except Exception as e:
                logger.warning(f"RCSB fetch failed ({url}): {e}")
                continue

        # WP-4: never fabricate a structure — fail loudly instead.
        raise SourceUnavailable(pdb_id, urls)

    @staticmethod
    def _acquire_afdb(candidate_id: str):
        """Fetch an AlphaFold prediction via metadata-first acquisition.

        P0 item 11: the hardcoded model_v4 file URLs are gone. The
        acquisition now queries /api/prediction/{uniprot}, selects the
        entry by `modelEntityId` (never data[0]) and `latestVersion`,
        downloads its `cifUrl`, and detects the format from the
        content. JSON is never accepted as a structure; every failure
        fails closed through SourceUnavailable with the AFDB record
        attached (including migrated identifiers, item 14).
        """
        try:
            pred = fetch_prediction(candidate_id)
        except AFDBUnavailable as exc:
            sources = [exc.record.get("api_url"),
                       exc.record.get("cif_url") or
                       "https://alphafold.ebi.ac.uk/files/"]
            # Attach the structured AFDB record (migrated identifiers,
            # fields present, versions) for the caller to log.
            exc.sources = [s for s in sources if s]
            raise SourceUnavailable(candidate_id, exc.sources) from exc
        fmt = "cif" if pred.fmt == "cif" else "pdb"
        logger.info(
            "AFDB Success: %s v%s (%s, %d bytes, sha256 %s...)",
            pred.model_entity_id, pred.latest_version, fmt,
            len(pred.content), pred.sha256[:12])
        return (pred.content, candidate_id, fmt)
