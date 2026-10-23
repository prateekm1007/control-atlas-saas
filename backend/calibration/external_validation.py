"""External engine validation harness skeleton.

Compares Toscanini bond/angle/Ramachandran/clash outlier calls against
MolProbity and wwPDB reference flags.

Roadmap requirements:
- >= 200 structures
- >= 95% agreement for each of 4 outlier categories
- Documented triage for every disagreement
- OD-16: no v25 cohort execution or threshold changes

This harness is a SKELETON. It defines the interface and data structures
but does not perform actual MolProbity comparison (requires external tools).
"""
from __future__ import annotations

import json
import hashlib
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class StructureValidation:
    """Validation result for a single structure."""
    pdb_id: str
    coordinate_sha256: str
    toscanini_bond_outliers: int
    toscanini_angle_outliers: int
    toscanini_rama_outliers: int
    toscanini_clash_outliers: int
    molprobity_bond_outliers: Optional[int] = None
    molprobity_angle_outliers: Optional[int] = None
    molprobity_rama_outliers: Optional[int] = None
    molprobity_clash_outliers: Optional[int] = None
    agreement_bond: Optional[bool] = None
    agreement_angle: Optional[bool] = None
    agreement_rama: Optional[bool] = None
    agreement_clash: Optional[bool] = None
    disagreement_triage: Optional[str] = None


@dataclass
class ValidationReport:
    """Aggregate validation report."""
    n_structures: int = 0
    bond_agreement_rate: float = 0.0
    angle_agreement_rate: float = 0.0
    rama_agreement_rate: float = 0.0
    clash_agreement_rate: float = 0.0
    gate_met: bool = False
    disagreements: list = field(default_factory=list)
    structures: list = field(default_factory=list)


CATEGORIES = ("bond", "angle", "rama", "clash")
MIN_STRUCTURES = 200
MIN_AGREEMENT = 0.95


def compute_agreement(toscanini: int, reference: int) -> bool:
    """Agreement: both agree on outlier count within tolerance."""
    return toscanini == reference


def build_report(validations: list[StructureValidation]) -> ValidationReport:
    report = ValidationReport()
    report.n_structures = len(validations)
    if not validations:
        return report

    for cat in CATEGORIES:
        agreements = 0
        for v in validations:
            attr = f"agreement_{cat}"
            if getattr(v, attr) is True:
                agreements += 1
        rate = agreements / len(validations)
        setattr(report, f"{cat}_agreement_rate", round(rate, 4))

    report.gate_met = (
        report.n_structures >= MIN_STRUCTURES
        and all(
            getattr(report, f"{cat}_agreement_rate") >= MIN_AGREEMENT
            for cat in CATEGORIES
        )
    )

    for v in validations:
        for cat in CATEGORIES:
            if getattr(v, f"agreement_{cat}") is False:
                report.disagreements.append({
                    "pdb_id": v.pdb_id,
                    "category": cat,
                    "toscanini": getattr(v, f"toscanini_{cat}_outliers"),
                    "reference": getattr(v, f"molprobity_{cat}_outliers"),
                    "triage": v.disagreement_triage,
                })

    report.structures = [
        {
            "pdb_id": v.pdb_id,
            "coordinate_sha256": v.coordinate_sha256,
            "toscanini": {
                "bond": v.toscanini_bond_outliers,
                "angle": v.toscanini_angle_outliers,
                "rama": v.toscanini_rama_outliers,
                "clash": v.toscanini_clash_outliers,
            },
            "molprobity": {
                "bond": v.molprobity_bond_outliers,
                "angle": v.molprobity_angle_outliers,
                "rama": v.molprobity_rama_outliers,
                "clash": v.molprobity_clash_outliers,
            },
        }
        for v in validations
    ]
    return report


def save_report(report: ValidationReport, path: Path) -> None:
    path.write_text(json.dumps({
        "n_structures": report.n_structures,
        "min_structures": MIN_STRUCTURES,
        "min_agreement": MIN_AGREEMENT,
        "bond_agreement_rate": report.bond_agreement_rate,
        "angle_agreement_rate": report.angle_agreement_rate,
        "rama_agreement_rate": report.rama_agreement_rate,
        "clash_agreement_rate": report.clash_agreement_rate,
        "gate_met": report.gate_met,
        "disagreement_count": len(report.disagreements),
        "disagreements": report.disagreements,
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    print("External validation harness skeleton.")
    print(f"Minimum structures: {MIN_STRUCTURES}")
    print(f"Minimum agreement: {MIN_AGREEMENT}")
    print("Status: NOT_EXECUTED (requires >= 200 structures + MolProbity)")
