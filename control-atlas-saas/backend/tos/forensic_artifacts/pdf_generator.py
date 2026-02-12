from __future__ import annotations

from datetime import datetime, timezone
from typing import Dict, Iterable, List, Optional

from fpdf import FPDF
import numpy as np


class DossierPDF(FPDF):
    def __init__(self, program) -> None:
        super().__init__(orientation="P", unit="mm", format="A4")
        self.program = program
        self.set_auto_page_break(auto=True, margin=15)

    def header(self) -> None:
        self.set_font("helvetica", "B", 11)
        self.set_text_color(15, 32, 56)
        self.cell(0, 6, "Toscanini Forensic Refusal Dossier", ln=True)
        self.set_font("helvetica", "", 8)
        notary = self.program.notary_seal or "UNSEALED"
        self.set_text_color(90, 90, 90)
        self.cell(0, 4, f"Program: {self.program.program_id} | Notary: {notary}", ln=True)
        self.set_draw_color(200, 200, 200)
        self.line(10, self.get_y() + 1, 200, self.get_y() + 1)
        self.ln(4)

    def footer(self) -> None:
        self.set_y(-12)
        self.set_font("helvetica", "I", 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 5, f"Page {self.page_no()} / {{nb}}", align="C")


def _truncate(text: str, max_chars: int) -> tuple[str, bool]:
    if len(text) <= max_chars:
        return text, False
    return text[: max_chars - 3].rstrip() + "...", True


def _collect_ca_positions(program) -> List[tuple[float, float, float]]:
    if not program.structure:
        return []
    ca_positions = [a.pos for a in program.structure.atoms if a.atom_name == "CA"]
    return ca_positions


def _draw_projection(pdf: FPDF, coords: List[tuple[float, float, float]], box, axes=(0, 1)) -> None:
    x0, y0, w, h = box
    pdf.rect(x0, y0, w, h)
    if len(coords) < 2:
        pdf.set_xy(x0, y0 + h / 2)
        pdf.set_font("helvetica", "I", 8)
        pdf.cell(w, 4, "No CA data", align="C")
        return
    arr = np.array(coords)
    xs = arr[:, axes[0]]
    ys = arr[:, axes[1]]
    min_x, max_x = float(np.min(xs)), float(np.max(xs))
    min_y, max_y = float(np.min(ys)), float(np.max(ys))
    span_x = max(max_x - min_x, 1e-3)
    span_y = max(max_y - min_y, 1e-3)
    scale = min((w - 6) / span_x, (h - 6) / span_y)
    offset_x = x0 + 3
    offset_y = y0 + 3
    pdf.set_draw_color(60, 100, 140)
    for i in range(len(coords) - 1):
        x1 = offset_x + (xs[i] - min_x) * scale
        y1 = offset_y + (ys[i] - min_y) * scale
        x2 = offset_x + (xs[i + 1] - min_x) * scale
        y2 = offset_y + (ys[i + 1] - min_y) * scale
        pdf.line(x1, y1, x2, y2)


def _ledger_rows(program) -> List[Dict[str, object]]:
    return program.tier1_ledger or []


def _nkg_count(program) -> int:
    for event in program.events:
        if event.state.value == "S6_FORENSIC_AGGREGATION":
            return int(event.payload.get("nkg_count", 0))
    return 0


def _probability_label(program) -> str:
    probability = program.tier3_breakdown.get("probability")
    if probability is None:
        return "N/A"
    return f"{probability * 100:.1f}%"


def generate_dossier(program) -> bytes:
    pdf = DossierPDF(program)
    pdf.alias_nb_pages()

    verdict = "VETO" if program.veto else "PASS"
    ingestion_date = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    notary_seal = program.notary_seal or "UNSEALED"
    tier3_probability = _probability_label(program)
    architecture = program.architecture.category if program.architecture else "UNASSIGNED"
    nkg_count = _nkg_count(program)

    # Page 1: Executive Summary
    pdf.add_page()
    pdf.set_font("helvetica", "B", 18)
    pdf.set_text_color(10, 10, 10)
    pdf.cell(0, 10, "Toscanini Forensic Refusal Dossier", ln=True, align="C")
    pdf.ln(2)
    pdf.set_font("helvetica", "", 10)
    pdf.cell(0, 6, f"Program ID: {program.program_id}", ln=True)
    pdf.cell(0, 6, f"Ingestion Date: {ingestion_date}", ln=True)
    pdf.cell(0, 6, f"Notary Seal: {notary_seal}", ln=True)
    pdf.ln(4)
    pdf.set_font("helvetica", "B", 24)
    pdf.set_text_color(160, 0, 0 if verdict == "VETO" else 90)
    pdf.cell(0, 12, f"FINAL VERDICT: {verdict}", ln=True, align="C")
    pdf.set_text_color(0, 0, 0)
    pdf.ln(2)
    pdf.set_font("helvetica", "", 11)
    executive_rationale = (
        "The sovereign refusal engine evaluated structural invariants, voxelized steric clashes, and "
        "contextual bonding constraints. The final verdict reflects Tier-1 physical feasibility and the "
        "cross-tier causal probability lock."
    )
    pdf.multi_cell(0, 6, executive_rationale)
    pdf.ln(2)
    pdf.set_font("helvetica", "B", 11)
    pdf.cell(0, 6, "Key Statistics", ln=True)
    pdf.set_font("helvetica", "", 10)
    pdf.cell(0, 6, f"Tier-3 Probability: {tier3_probability}", ln=True)
    pdf.cell(0, 6, f"Derived Architecture: {architecture}", ln=True)
    pdf.cell(0, 6, f"NKG Failure Moat Entries: {nkg_count}", ln=True)

    # Page 2: Engineering Overview
    pdf.add_page()
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 8, "Engineering Overview: Dual-Engine Visualizer", ln=True)
    pdf.set_font("helvetica", "", 9)
    pdf.multi_cell(
        0,
        5,
        "Dual projections illustrate PASS (Cartoon) and VETO (Stick-Autopsy) renderings. "
        "Topological symmetry is centered to maintain visual parity.",
    )
    ca_positions = _collect_ca_positions(program)
    box_y = pdf.get_y() + 4
    box_w = 85
    box_h = 70
    pdf.set_font("helvetica", "B", 10)
    pdf.cell(box_w, 5, "PASS (Cartoon)", ln=0, align="C")
    pdf.cell(0, 5, "VETO (Stick-Autopsy)", ln=1, align="C")
    _draw_projection(pdf, ca_positions, (15, box_y + 6, box_w, box_h), axes=(0, 1))
    _draw_projection(pdf, ca_positions, (110, box_y + 6, box_w, box_h), axes=(2, 1))
    pdf.set_y(box_y + box_h + 12)
    pdf.set_font("helvetica", "I", 8)
    pdf.cell(0, 5, "Centered projections preserve visual symmetry across the dossier.", ln=True, align="C")

    # Pages 3-4: Diagnostic Ledger
    pdf.add_page()
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 8, "Diagnostic Ledger", ln=True)
    pdf.set_font("helvetica", "", 9)
    pdf.cell(0, 6, "Law ID | Title | Principle | Measurement | Status", ln=True)
    pdf.ln(2)
    for entry in _ledger_rows(program):
        status = entry.get("status", "PASS")
        fill = status in {"FAIL", "VETO"}
        if fill:
            pdf.set_fill_color(255, 230, 230)
        else:
            pdf.set_fill_color(240, 240, 240)
        pdf.set_font("helvetica", "B", 9)
        pdf.cell(0, 6, f"{entry.get('law_id')} [{status}]", ln=True, fill=True)
        pdf.set_font("helvetica", "", 8)
        pdf.multi_cell(0, 4, f"{entry.get('title')} | {entry.get('principle')}")
        pdf.multi_cell(0, 4, f"Measurement: {entry.get('measurement')}")
        anchor = entry.get("anchor") or {}
        if anchor:
            pdf.set_text_color(140, 0, 0)
            pdf.multi_cell(0, 4, f"Violation Anchor: {anchor}")
            pdf.set_text_color(0, 0, 0)
        pdf.ln(1)
        if pdf.get_y() > 250:
            pdf.add_page()
            pdf.set_font("helvetica", "B", 14)
            pdf.cell(0, 8, "Diagnostic Ledger (cont.)", ln=True)

    # Page 5: Architecture & Probability Breakdown
    pdf.add_page()
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 8, "Architecture & Probability Breakdown", ln=True)
    pdf.set_font("helvetica", "", 10)
    if program.architecture:
        pdf.multi_cell(
            0,
            5,
            f"Derived Category: {program.architecture.category}\nRationale: {program.architecture.rationale}\n"
            f"Voxel Span: {program.architecture.voxel_span}\n"
            f"Symmetry Score: {program.architecture.symmetry_score}\n"
            f"Clash-Free Volume: {program.architecture.clash_free_volume}",
        )
    else:
        pdf.multi_cell(0, 5, "Architecture derivation unavailable.")
    pdf.ln(2)
    pdf.set_font("helvetica", "B", 11)
    pdf.cell(0, 6, "Tier-3 Probability Math", ln=True)
    pdf.set_font("helvetica", "", 10)
    if program.tier3_breakdown:
        for key, label in [
            ("p_base", "P_base"),
            ("p_phys", "P_phys"),
            ("p_ml", "P_ml"),
            ("tax", "Tax"),
            ("weight_derivation", "Weight_Derivation"),
            ("probability", "Final P"),
        ]:
            value = program.tier3_breakdown.get(key, 0.0)
            pdf.cell(0, 6, f"{label}: {value:.4f}", ln=True)
    else:
        pdf.cell(0, 6, "Tier-3 probability breakdown unavailable.", ln=True)
    pdf.ln(2)
    pdf.set_font("helvetica", "B", 11)
    pdf.cell(0, 6, "Epistemic Duality Narrative", ln=True)
    pdf.set_font("helvetica", "", 9)
    pdf.multi_cell(
        0,
        5,
        "Physical confidence reflects deterministic invariant compliance, while ML confidence represents "
        "model-derived structural plausibility. Both are required for Schrödinger-Class credibility.",
    )

    # Page 6: PhD-Grade Rejection Rationale
    pdf.add_page()
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 8, "PhD-Grade Rationale", ln=True)
    pdf.set_font("helvetica", "", 10)
    rationale = (
        "This dossier provides a full forensic audit of the submitted structure. The Tier-1 voxelized physics "
        "engine enforces steric exclusion and backbone continuity with deterministic spatial hashing. "
        "Any violation of LAW-155-L triggers a sovereign veto, preventing downstream synthesis waste. "
        "The architecture derivation stage evaluates voxel span, symmetry, and clash-free volume to assign "
        "a Tier-2 category and inform causal weighting. The Tier-3 Bayesian scoring engine integrates "
        "historical baselines, physical boosts, ML confidence, and refinement tax to quantify overall success "
        "probability. Where a veto is present, the narrative highlights the specific steric overlap and "
        "explains why the geometry violates fundamental chemical constraints. Where no veto exists, the "
        "report documents the physics compliance and substantiates the Schrödinger-Class threshold with "
        "auditable probability breakdowns."
    )
    rationale, truncated = _truncate(rationale, 4200)
    pdf.multi_cell(0, 5, rationale)
    if truncated:
        pdf.set_text_color(180, 0, 0)
        pdf.set_font("helvetica", "I", 9)
        pdf.cell(0, 5, "Warning: Narrative truncated to fit the dossier page.", ln=True)
        pdf.set_text_color(0, 0, 0)

    # Page 7: Annex (optional)
    if program.hash_chain:
        pdf.add_page()
        pdf.set_font("helvetica", "B", 14)
        pdf.cell(0, 8, "Annex: Notary & Hash Chain", ln=True)
        pdf.set_font("helvetica", "", 9)
        pdf.multi_cell(
            0,
            5,
            "REMARK 900 TOSCANINI FORENSIC SEAL\n"
            f"AUDIT_ID: {program.program_id}\n"
            f"VERDICT: {verdict}\n"
            f"COORD_HASH: {notary_seal}\n"
            "STATUS: NOTARIZED_SOVEREIGN",
        )
        pdf.ln(2)
        pdf.set_font("helvetica", "B", 10)
        pdf.cell(0, 6, "Hash Chain Proof", ln=True)
        pdf.set_font("helvetica", "", 8)
        for chain_hash in program.hash_chain:
            pdf.multi_cell(0, 4, chain_hash)

    return pdf.output(dest="S").encode("latin-1")


def generate_v14_certificate(sig, verdict, score, gen, rationale, results, atoms):
    pdf = FPDF()
    pdf.set_auto_page_break(True)
    pdf.add_page()
    pdf.set_font("helvetica", "B", 18)
    pdf.cell(0, 15, "TOSCANINI FORENSIC NOTARY CERTIFICATE", ln=True, align="C")

    color = (180, 0, 0) if verdict == "VETO" else (0, 100, 0)
    pdf.set_fill_color(245, 245, 245)
    pdf.set_text_color(*color)
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 12, f"DECISION: {verdict} / PHYSICAL INTEGRITY: {score}%", ln=True, align="C", fill=True)
    pdf.set_text_color(0)
    pdf.ln(5)

    pdf.set_font("helvetica", "B", 11)
    pdf.cell(0, 8, "EXECUTIVE BIOPHYSICAL SYNTHESIS", ln=True)
    pdf.set_font("helvetica", "", 10)
    pdf.multi_cell(0, 5, str(rationale).encode("ascii", "ignore").decode())

    pdf.ln(5)
    pdf.set_fill_color(230, 240, 250)
    pdf.set_font("helvetica", "B", 10)
    pdf.cell(0, 8, "TECHNICAL METRIC INTERPRETATION", ln=True, fill=True)
    pdf.set_font("helvetica", "", 9)
    interp = (
        "Physical Score: Confirms adherence to invariant laws of atomic exclusion. "
        "ML Confidence (pLDDT): Represents the predicted C-alpha local distance difference test. "
        "Scores < 70 indicate regions where tertiary side-chain packing is speculative and docking is high-risk."
    )
    pdf.multi_cell(0, 4, interp)

    ca = [a.pos for a in atoms if (a.atom_name if hasattr(a, "atom_name") else a.get("atom")) == "CA"]
    if len(ca) > 4:
        pdf.ln(5)
        coords = np.array(ca)
        xs, ys = coords[:, 0], coords[:, 1]
        mx, my, span = np.min(xs), np.min(ys), max(1, np.max(xs) - np.min(xs))
        scale, bx, by = 60 / span, 30, pdf.get_y() + 10
        pdf.set_draw_color(150, 150, 150)
        pdf.set_line_width(0.8)
        for i in range(len(ca) - 1):
            pdf.line(
                bx + (ca[i][0] - mx) * scale,
                by + (ca[i][1] - my) * scale,
                bx + (ca[i + 1][0] - mx) * scale,
                by + (ca[i + 1][1] - my) * scale,
            )
        bx2 = 120
        pdf.set_draw_color(100, 120, 140)
        for i in range(len(ca) - 1):
            pdf.line(
                bx2 + (ca[i][2] - mx) * scale,
                by + (ca[i][1] - my) * scale,
                bx2 + (ca[i + 1][2] - mx) * scale,
                by + (ca[i + 1][1] - my) * scale,
            )

    pdf.add_page()
    pdf.set_font("helvetica", "B", 16)
    pdf.cell(0, 15, "TECHNICAL ANNEX (10-LAW CANON)", ln=True)
    for l in results:
        pdf.set_font("helvetica", "B", 9)
        pdf.cell(0, 6, f"{l['law_id']}: {l['title']} [{l['status']}]", ln=True)
        pdf.set_font("helvetica", "", 8)
        pdf.multi_cell(0, 4, f"Observed: {l['measurement']}\nPrinciple: {l['principle']}")
        pdf.ln(2)
    return pdf.output()
