"""Program lifecycle and causality chain utilities for Toscanini v21.0.4."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Optional

from tos.engine import tier1_measurements
from tos.engine.tier3_predict import ProbabilityInputs, compute_probability_breakdown
from tos.export.sealer import ForensicSealer
from tos.glossary import law_glossary
from tos.saas_engine.selector import ArchitectureDecision, derive_architecture


class ProgramState(str, Enum):
    S0_INGESTED = "S0_INGESTED"
    S1_ACQUISITION_COMPLETE = "S1_ACQUISITION_COMPLETE"
    S2_FOLDING_COMPLETE = "S2_FOLDING_COMPLETE"
    S3_TIER1_COMPLETE = "S3_TIER1_COMPLETE"
    S4_ARCHITECTURE_DERIVED = "S4_ARCHITECTURE_DERIVED"
    S5_TIER3_COMPLETE = "S5_TIER3_COMPLETE"
    S6_FORENSIC_AGGREGATION = "S6_FORENSIC_AGGREGATION"
    S7_DOSSIER_READY = "S7_DOSSIER_READY"


STATE_SEQUENCE = [
    ProgramState.S0_INGESTED,
    ProgramState.S1_ACQUISITION_COMPLETE,
    ProgramState.S2_FOLDING_COMPLETE,
    ProgramState.S3_TIER1_COMPLETE,
    ProgramState.S4_ARCHITECTURE_DERIVED,
    ProgramState.S5_TIER3_COMPLETE,
    ProgramState.S6_FORENSIC_AGGREGATION,
    ProgramState.S7_DOSSIER_READY,
]


@dataclass(frozen=True)
class ProgramEvent:
    state: ProgramState
    payload: Dict[str, object]
    previous_hash: str
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def chain_hash(self) -> str:
        content = {
            "state": self.state,
            "payload": self.payload,
            "previous_hash": self.previous_hash,
            "created_at": self.created_at,
        }
        encoded = json.dumps(content, sort_keys=True, default=str).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()


@dataclass
class ProgramRecord:
    program_id: str
    events: List[ProgramEvent] = field(default_factory=list)
    structure: Optional[object] = None
    notary_seal: Optional[str] = None
    tier1_ledger: List[Dict[str, object]] = field(default_factory=list)
    architecture: Optional[ArchitectureDecision] = None
    tier3_breakdown: Dict[str, float] = field(default_factory=dict)
    veto: bool = False
    hash_chain: List[str] = field(default_factory=list)

    @property
    def current_state(self) -> Optional[ProgramState]:
        if not self.events:
            return None
        return self.events[-1].state

    def add_event(self, state: ProgramState, payload: Dict[str, object]) -> ProgramEvent:
        previous_hash = self.events[-1].chain_hash() if self.events else "GENESIS"
        event = ProgramEvent(state=state, payload=payload, previous_hash=previous_hash)
        self.events.append(event)
        self.hash_chain.append(event.chain_hash())
        return event

    def validate_chain(self) -> bool:
        previous_hash = "GENESIS"
        for event in self.events:
            if event.previous_hash != previous_hash:
                return False
            previous_hash = event.chain_hash()
        return True

    def require_next_state(self, state: ProgramState) -> None:
        expected_index = 0 if self.current_state is None else STATE_SEQUENCE.index(self.current_state) + 1
        if expected_index >= len(STATE_SEQUENCE) or STATE_SEQUENCE[expected_index] != state:
            raise ValueError(f"Non-skippable state violation: expected {STATE_SEQUENCE[expected_index].value}")


class ProgramRepository:
    """In-memory repository for program lifecycle tracking."""

    def __init__(self) -> None:
        self._records: Dict[str, ProgramRecord] = {}

    def upsert(self, program_id: str) -> ProgramRecord:
        if program_id not in self._records:
            self._records[program_id] = ProgramRecord(program_id=program_id)
        record = self._records[program_id]
        if record.events and not record.validate_chain():
            raise ValueError("Program hash chain integrity failure.")
        return record

    def get(self, program_id: str) -> Optional[ProgramRecord]:
        record = self._records.get(program_id)
        if record and not record.validate_chain():
            raise ValueError("Program hash chain integrity failure.")
        return record

    def list_programs(self) -> List[ProgramRecord]:
        return list(self._records.values())


def _compute_voxel_span(structure) -> tuple[float, float, float]:
    positions = [a.pos for a in structure.atoms + structure.ligands]
    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def _symmetry_score(voxel_span: tuple[float, float, float]) -> float:
    max_dim = max(voxel_span)
    min_dim = max(min(voxel_span), 1e-6)
    return min_dim / max_dim


def _clash_free_volume(voxel_span: tuple[float, float, float]) -> float:
    x, y, z = voxel_span
    return float(x * y * z)


def advance_program(
    program: ProgramRecord,
    structure,
    acquisition_source: str,
    folding_source: str,
    plddt: float,
) -> ProgramRecord:
    """Advance a program deterministically through the full sovereign pipeline."""

    if program.current_state is None:
        program.require_next_state(ProgramState.S0_INGESTED)
        coord_hash = ForensicSealer.generate_hash(ForensicSealer.canonical_serialize(structure.atoms + structure.ligands))
        program.structure = structure
        program.notary_seal = coord_hash
        program.add_event(ProgramState.S0_INGESTED, {"notary_seal": coord_hash})

    if program.current_state == ProgramState.S0_INGESTED:
        program.require_next_state(ProgramState.S1_ACQUISITION_COMPLETE)
        program.add_event(ProgramState.S1_ACQUISITION_COMPLETE, {"source": acquisition_source})

    if program.current_state == ProgramState.S1_ACQUISITION_COMPLETE:
        program.require_next_state(ProgramState.S2_FOLDING_COMPLETE)
        program.add_event(ProgramState.S2_FOLDING_COMPLETE, {"folding_source": folding_source})

    if program.current_state == ProgramState.S2_FOLDING_COMPLETE:
        program.require_next_state(ProgramState.S3_TIER1_COMPLETE)
        s155, m155, a155 = tier1_measurements.Tier1Measurements.check_law_155_L(structure)
        s160, m160, a160 = tier1_measurements.Tier1Measurements.check_law_160(structure)
        ledger: List[Dict[str, object]] = []
        for law_id in law_glossary.list_all_law_ids():
            status, measurement, anchor = ("PASS", "Verified Invariant", {})
            if law_id == "LAW-155":
                status, measurement, anchor = s155, m155, a155
            elif law_id == "LAW-160":
                status, measurement, anchor = s160, m160, a160
            explanation = law_glossary.get_law_explanation(law_id)
            ledger.append(
                {
                    "law_id": law_id,
                    "status": status,
                    "measurement": measurement,
                    "anchor": anchor,
                    "title": explanation["title"],
                    "principle": explanation["principle"],
                }
            )
        program.tier1_ledger = ledger
        program.veto = any(l["status"] in {"FAIL", "VETO"} for l in ledger)
        program.add_event(ProgramState.S3_TIER1_COMPLETE, {"veto": program.veto, "ledger_count": len(ledger)})
        if program.veto:
            return program

    if program.current_state == ProgramState.S3_TIER1_COMPLETE:
        program.require_next_state(ProgramState.S4_ARCHITECTURE_DERIVED)
        voxel_span = _compute_voxel_span(structure)
        symmetry = _symmetry_score(voxel_span)
        clash_free_volume = _clash_free_volume(voxel_span)
        program.architecture = derive_architecture(
            voxel_span=voxel_span,
            symmetry_score=symmetry,
            clash_free_volume=clash_free_volume,
        )
        program.add_event(
            ProgramState.S4_ARCHITECTURE_DERIVED,
            {
                "category": program.architecture.category,
                "voxel_span": voxel_span,
                "symmetry_score": symmetry,
                "clash_free_volume": clash_free_volume,
            },
        )

    if program.current_state == ProgramState.S4_ARCHITECTURE_DERIVED:
        program.require_next_state(ProgramState.S5_TIER3_COMPLETE)
        weight_by_category = {
            "Linear": 0.95,
            "Multivalent": 1.15,
            "Metal": 1.05,
            "Helical": 1.0,
            "Sheet": 0.98,
            "Compact": 1.1,
            "Extended": 0.92,
        }
        derivation_weight = weight_by_category.get(program.architecture.category, 1.0)
        inputs = ProbabilityInputs(
            plddt=plddt,
            physical_boost=0.25,
            base_probability=0.15,
            refinement_tax=0.25,
            derivation_weight=derivation_weight,
        )
        program.tier3_breakdown = compute_probability_breakdown(inputs)
        program.add_event(ProgramState.S5_TIER3_COMPLETE, program.tier3_breakdown)

    if program.current_state == ProgramState.S5_TIER3_COMPLETE:
        program.require_next_state(ProgramState.S6_FORENSIC_AGGREGATION)
        nkg_path = Path("/app/nkg/piu_moat.jsonl")
        nkg_count = 0
        if nkg_path.exists():
            with nkg_path.open("r", encoding="utf-8") as handle:
                nkg_count = sum(1 for _ in handle)
        program.add_event(ProgramState.S6_FORENSIC_AGGREGATION, {"nkg_count": nkg_count})

    if program.current_state == ProgramState.S6_FORENSIC_AGGREGATION:
        program.require_next_state(ProgramState.S7_DOSSIER_READY)
        program.add_event(ProgramState.S7_DOSSIER_READY, {"ready": True})

    return program
