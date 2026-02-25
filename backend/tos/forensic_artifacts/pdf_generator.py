from fpdf import FPDF
from .confidence_viz import generate_confidence_bar
from .structure_viz import generate_structure_render
from .residue_heatmap import generate_residue_heatmap_png
from datetime import datetime
from .pdf_theme import PDF_COLORS, PDF_LAYOUT
from ..governance.constants import LAW_105_THRESHOLD
from ..utils.type_guards import force_bytes
from ..governance.station_sop import (
    STATION_METADATA, BAYESIAN_FORMULA, LAW_CANON_HASH
)
from ..governance.modality_matrix import compute_matrix_hash


class ToscaniniDossier(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 8)
        self.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
        self.cell(w=0, h=10, txt=f"TOSCANINI STRUCTURAL GOVERNANCE DOSSIER v{STATION_METADATA['version']}", border=0, align="L")
        self.ln(10)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 7)
        self.set_text_color(*PDF_COLORS["TEXT_MUTED"])
        self.cell(0, 10, f"Engine: v{STATION_METADATA['version']} // Canon: {LAW_CANON_HASH[:12]}... // Matrix: {compute_matrix_hash()[:12]}... // Page {self.page_no()}", align="R")

    def section_title(self, text):
        self.set_font("Helvetica", "B", 11)
        self.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        self.cell(0, 8, text.upper(), ln=True)
        self.set_draw_color(0, 0, 0)
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(4)

    def verdict_badge(self, verdict: str, score: int, coverage: float):
        """Large centered verdict badge for executive cognition."""
        color_map = {
            "PASS": PDF_COLORS["PASS"],
            "VETO": PDF_COLORS["VETO"],
            "INDETERMINATE": PDF_COLORS["INDETERMINATE"],
        }
        badge_color = color_map.get(verdict, PDF_COLORS["INDETERMINATE"])
        badge_width = 190

        # Badge background
        self.set_fill_color(*badge_color)
        self.set_text_color(*PDF_COLORS["TEXT_WHITE"])
        self.set_font("Helvetica", "B", 32)
        self.cell(badge_width, 26, verdict, border=0, ln=True, align="C", fill=True)

        # Score line below badge
        self.set_fill_color(*PDF_COLORS["BG_LIGHT"])
        self.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        self.set_font("Helvetica", "B", 11)
        self.cell(badge_width, 10,
                  f"Deterministic Score: {score}/100   |   Coverage: {coverage}%",
                  border=0, ln=True, align="C", fill=True)
        self.ln(6)

    def safe(self, text):
        return str(text).encode('latin-1', 'replace').decode('latin-1')

    def red_banner(self, text):
        """Regulator Finding #1: High-contrast epistemic warning."""
        self.set_fill_color(*PDF_COLORS["BANNER_VETO"])
        self.set_text_color(*PDF_COLORS["TEXT_WHITE"])
        self.set_font("Helvetica", "B", 10)
        self.cell(190, 10, self.safe(f"  {text}"), border=1, ln=True, fill=True, align="L")
        self.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        self.ln(3)

    def amber_banner(self, text):
        """Amber warning for advisory notices."""
        self.set_fill_color(*PDF_COLORS["BANNER_ALERT"])
        self.set_text_color(*PDF_COLORS["TEXT_WHITE"])
        self.set_font("Helvetica", "B", 9)
        self.cell(190, 9, self.safe(f"  {text}"), border=1, ln=True, fill=True, align="L")
        self.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        self.ln(3)

    def section_divider(self, text, r, g, b):
        """Regulator Finding #4: Visually distinct section headers."""
        self.set_fill_color(r, g, b)
        self.set_text_color(*PDF_COLORS["TEXT_WHITE"])
        self.set_font("Helvetica", "B", 10)
        self.cell(190, 9, self.safe(f"  {text}"), border=0, ln=True, fill=True, align="L")
        self.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        self.ln(3)


def _extract_ca_atoms(pdb_b64):
    """Extract CA atom coordinates and B-factors from base64-encoded PDB.
    Used for 3D structure visualization in the forensic dossier.
    Returns: (list of (x,y,z) tuples, list of confidence floats)"""
    import base64 as _b64
    coords, confidences = [], []
    try:
        pdb_text = _b64.b64decode(pdb_b64).decode("utf-8", errors="ignore")
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and len(line) >= 54:
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                    try:
                        bfactor = float(line[60:66])
                        confidences.append(bfactor)
                    except (ValueError, IndexError):
                        confidences.append(50.0)
    except Exception:
        pass
    return coords, confidences



def _find_low_regions(confidences, threshold=50.0, min_length=3):
    """Find contiguous stretches of residues below threshold. Returns 1-indexed tuples."""
    regions, in_r, start = [], False, 0
    for i, s in enumerate(confidences):
        if s < threshold:
            if not in_r: start, in_r = i, True
        elif in_r:
            if (i - start) >= min_length: regions.append((start + 1, i))
            in_r = False
    if in_r and (len(confidences) - start) >= min_length:
        regions.append((start + 1, len(confidences)))
    return regions

def generate_v21_dossier(payload):
    try:
        pdf = ToscaniniDossier()
        pdf.set_auto_page_break(True, margin=15)
        pdf.set_compression(True)
        w = 190
        v = payload.get('verdict', {})
        ss = payload.get("characterization", {})
        prov = payload.get("provenance", {})
        conf_meta = payload.get("confidence_meta", {})
        gov = payload.get("governance", {})
        
        # PDF Metadata for searchability and institutional indexing
        audit_id = gov.get("audit_id", "UNKNOWN")
        pdf.set_title(f"Toscanini Structural Governance Certification - {audit_id}")
        pdf.set_author("Toscanini Structural Governance Engine")
        pdf.set_subject("Deterministic Structural Admissibility Report")
        pdf.set_keywords("protein structure, governance, certification, forensic, AlphaFold")
        math_data = payload.get("strategic_math", {})
        laws = payload.get('tier1', {}).get('laws', [])

        is_experimental = ss.get("source_type") != "predicted"
        verdict_binary = v.get("binary", "ERROR")
        coverage = v.get("coverage_pct", 0)

        # ══════ PAGE 1: EXECUTIVE DETERMINATION ══════
        pdf.add_page()
        pdf.section_title("Executive Determination")

        # Source Attribution Header
        source_label = f"EXPERIMENTAL ({ss.get('source_type', 'N/A').upper()})" if is_experimental else "ML-PREDICTED (ALPHAFOLD)"
        pdf.set_font("Helvetica", "B", 9)
        pdf.set_fill_color(*PDF_COLORS["TEXT_SECONDARY"])
        pdf.set_text_color(*PDF_COLORS["TEXT_WHITE"])
        pdf.cell(w, 7, f"DATA REGIME: {source_label}", border=0, ln=True, fill=True, align="C")
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        pdf.ln(4)

        # VERDICT BADGE - Large centered institutional panel
        det_score = v.get("deterministic_score", 0)
        pdf.verdict_badge(verdict_binary, det_score, coverage)

        # METRIC CARDS - Institutional cards with left accent strip
        card_w = 56
        card_gap = 7
        card_total_h = 27
        start_x = PDF_LAYOUT["MARGIN_LEFT"]

        det_pass = v.get('det_passed', 0)
        det_total_v = v.get('det_total', 12)
        cov_val = coverage
        pri_val = payload.get('tier3', {}).get('probability', 0)

        if cov_val >= 70:
            cov_sub, cov_col = "SUFFICIENT", PDF_COLORS["PASS"]
        elif cov_val >= 50:
            cov_sub, cov_col = "MARGINAL", PDF_COLORS["ALERT"]
        else:
            cov_sub, cov_col = "INSUFFICIENT", PDF_COLORS["INDETERMINATE"]

        if det_pass == det_total_v:
            comp_sub, comp_col = "ALL LAWS MET", PDF_COLORS["PASS"]
        elif det_pass >= det_total_v * 0.8:
            comp_sub, comp_col = "MINOR FLAGS", PDF_COLORS["ALERT"]
        else:
            comp_sub, comp_col = "VIOLATIONS", PDF_COLORS["VETO"]

        if pri_val >= 70:
            pri_sub, pri_col = "HIGH PRIORITY", PDF_COLORS["PASS"]
        elif pri_val >= 40:
            pri_sub, pri_col = "MODERATE", PDF_COLORS["ALERT"]
        else:
            pri_sub, pri_col = "LOW / GATED", PDF_COLORS["INDETERMINATE"]

        cards = [
            ("DETERMINISTIC COMPLIANCE", "%d / %d" % (det_pass, det_total_v), comp_sub, comp_col),
            ("RELIABILITY COVERAGE",     "%.1f%%" % cov_val,                  cov_sub,  cov_col),
            ("PRIORITIZATION INDEX",     "%s%%" % str(pri_val),               pri_sub,  pri_col),
        ]

        card_y = pdf.get_y()

        for idx, (label, value, subtext, val_color) in enumerate(cards):
            cx = start_x + idx * (card_w + card_gap)
            accent_w = 2.5

            # Left accent strip (colored indicator)
            pdf.set_fill_color(*val_color)
            pdf.rect(cx, card_y, accent_w, card_total_h, style='F')

            # Card body (single rect, no internal borders)
            body_x = cx + accent_w
            body_w = card_w - accent_w
            pdf.set_draw_color(215, 215, 215)
            pdf.set_line_width(0.25)
            pdf.set_fill_color(*PDF_COLORS["BG_LIGHTER"])
            pdf.rect(body_x, card_y, body_w, card_total_h, style='DF')

            # Label (top zone - grey text, 7pt)
            pdf.set_xy(body_x, card_y + 2)
            pdf.set_font("Helvetica", "B", 7)
            pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
            pdf.cell(body_w, 5, label, border=0, align="C")

            # Value (center zone - black, 18pt bold)
            pdf.set_xy(body_x, card_y + 8)
            pdf.set_font("Helvetica", "B", 18)
            pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
            pdf.cell(body_w, 12, value, border=0, align="C")

            # Subtext (bottom zone - colored status)
            pdf.set_xy(body_x, card_y + 20)
            pdf.set_font("Helvetica", "B", 8)
            pdf.set_text_color(*val_color)
            pdf.cell(body_w, 5, subtext, border=0, align="C")

        # Reset drawing state below cards
        pdf.set_draw_color(0, 0, 0)
        pdf.set_line_width(0.2)
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        pdf.set_y(card_y + card_total_h + 7)
        pdf.set_x(PDF_LAYOUT["MARGIN_LEFT"])

        # FULL-WIDTH pLDDT CONFIDENCE BAR
        try:
            bar_png = generate_confidence_bar(
                coverage, v.get('det_passed', 0), v.get('det_total', 12)
            )
            import tempfile as _tf_bar
            import os as _os_bar
            with _tf_bar.NamedTemporaryFile(suffix='.png', delete=False) as _tmp_bar:
                _tmp_bar.write(bar_png)
                _bar_path = _tmp_bar.name
            pdf.image(_bar_path, x=PDF_LAYOUT["MARGIN_LEFT"], w=190, h=30)
            _os_bar.unlink(_bar_path)
            pdf.ln(5)
        except Exception:
            pdf.set_font("Helvetica", "I", 8)
            pdf.cell(190, 6, "[Confidence distribution unavailable]",
                     ln=True, align="C")
            pdf.ln(2)

        # Obligation statement
        pdf.set_font("Helvetica", "B", 9)
        pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
        pdf.cell(w, 6,
                 "Independent experimental confirmation required prior to reliance on this certification.",
                 border=0, ln=True, align="C")
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        pdf.ln(4)

        # ═══ REGULATOR FIX #1: Coverage failure high-contrast warning ═══
        if coverage < LAW_105_THRESHOLD:
            pdf.red_banner(f"ADJUDICATION SUSPENDED -- Evidence Scope Below Threshold (LAW-105: {coverage}%)")
            pdf.set_font("Helvetica", "I", 8)
            pdf.multi_cell(w, 4, "The engine cannot certify structural compliance when reliability "
                           "coverage falls below the 70% epistemic threshold. This is not a declaration "
                           "of failure -- it reflects insufficient high-confidence coordinate data to "
                           "support deterministic governance.")
            pdf.ln(3)

        if verdict_binary == "VETO":
            det_fails = [l for l in laws if l.get('method') == 'deterministic' and l['status'] != 'PASS']
            fail_ids = ', '.join([l['law_id'] for l in det_fails])
            pdf.red_banner(f"DETERMINISTIC VETO: {fail_ids}")

        pdf.ln(3)
        # Interpretation Clause
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "INTERPRETATION OF VERDICT CLASSIFICATION", ln=True)
        pdf.set_font("Helvetica", "", 8)

        pdf.set_font("Helvetica", "B", 8); pdf.write(5, "PASS: "); pdf.set_font("Helvetica", "", 8)
        pdf.write(5, f"Physically compliant. All deterministic invariants satisfied. Coverage >= {LAW_105_THRESHOLD}%.\n")
        pdf.set_font("Helvetica", "B", 8); pdf.write(5, "VETO: "); pdf.set_font("Helvetica", "", 8)
        pdf.write(5, "Physically incompatible. Non-negotiable invariant violation detected.\n")
        pdf.set_font("Helvetica", "B", 8); pdf.write(5, "INDETERMINATE: "); pdf.set_font("Helvetica", "", 8)
        pdf.write(5, "Epistemic insufficiency. Coverage < 70% or low-confidence data. Not a failure declaration.\n")

        # ═══ REGULATOR FIX #2: Modality jurisdiction statement ═══
        pdf.ln(4)
        pdf.amber_banner("PIL-CAL-03: Modality-Aware Enforcement Active")
        pdf.set_font("Helvetica", "I", 8)
        pdf.multi_cell(w, 4, "Deterministic enforcement is conditioned by acquisition modality per PIL-CAL-03. "
                       "X-ray: LAW-100 advisory. Cryo-EM: LAW-100, LAW-170 advisory. "
                       "NMR: LAW-100, LAW-125, LAW-170 advisory. "
                       "Predicted (AlphaFold): all deterministic laws enforced at full jurisdiction.")

        # ══════ PAGE 2: SECTION I - DETERMINISTIC GOVERNANCE ══════
        pdf.add_page()

        # ═══ REGULATOR FIX #4: Distinct section header ═══
        pdf.section_divider("SECTION I: DETERMINISTIC GOVERNANCE (SCORE-DRIVING)", 20, 60, 120)

        pdf.set_font("Helvetica", "I", 8)
        pdf.multi_cell(w, 4, "Physical Refusal Gate. Failure of any deterministic law triggers VETO for predicted models. "
                       "For experimental structures, modality-specific reclassifications apply per PIL-CAL-03.")
        pdf.ln(2)
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_text_color(*PDF_COLORS["ALERT"])
        # CANONICAL ADJUDICATION HIERARCHY - DO NOT DUPLICATE ANYWHERE
        pdf.multi_cell(w, 4, f"Adjudication Hierarchy: Coverage gate (LAW-105 >= {LAW_105_THRESHOLD}%) precedes deterministic incompatibility checks. If coverage is insufficient, verdict is INDETERMINATE regardless of individual law compliance. This reflects epistemic insufficiency, not structural failure.")
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        pdf.ln(2)
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_text_color(*PDF_COLORS["ALERT"])

        det_laws = [l for l in laws if l.get('method') == 'deterministic']
        violations = [l for l in det_laws if l['status'] != 'PASS']
        passes = [l for l in det_laws if l['status'] == 'PASS']

        headers = ["ID", "Metric", "Class", "Observed", "Threshold", "Status"]
        widths = [15, 50, 30, 35, 35, 25]

        if violations:
            pdf.set_font("Helvetica", "B", 9)
            pdf.set_text_color(*PDF_COLORS["VETO"])
            pdf.cell(w, 8, f"DETECTION: {len(violations)} DETERMINISTIC VIOLATION(S)", ln=True)
            pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
            pdf.set_font("Helvetica", "B", 8)
            pdf.set_fill_color(*PDF_COLORS["BG_VETO_TABLE"])
            for i, h in enumerate(headers):
                pdf.cell(widths[i], 8, h, 1, 0, "C", True)
            pdf.ln(8)
            pdf.set_font("Helvetica", "", 7)
            for l in violations:
                for i, field in enumerate([l['law_id'], l['title'], "Deterministic",
                                           str(l['observed']), f"{l['operator']} {l['threshold']}", l['status']]):
                    pdf.cell(widths[i], 7, pdf.safe(field), 1, 0, "C" if i != 1 else "L")
                pdf.ln(7)
            pdf.ln(5)

        pdf.set_font("Helvetica", "B", 9)
        pdf.set_font("Helvetica", "B", 9)
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        pdf.cell(w, 8, "DETERMINISTIC COMPLIANCE LEDGER", ln=True)
        pdf.set_font("Helvetica", "B", 7)
        pdf.set_fill_color(*PDF_COLORS["BG_PASS_TABLE"])
        for i, h in enumerate(headers):
            pdf.cell(widths[i], 7, h, 1, 0, "C", True)
        pdf.ln(7)
        pdf.set_font("Helvetica", "", 7)
        for law in passes:
            row_fields = [
                law['law_id'],
                law['title'],
                "Deterministic",
                str(law.get('observed', 'N/A')),
                "%s %s" % (law.get('operator', ''), law.get('threshold', 'N/A')),
                "[PASS]"
            ]
            pdf.set_fill_color(*PDF_COLORS["BG_PASS_TABLE"])
            for i, field in enumerate(row_fields):
                pdf.cell(widths[i], 6, pdf.safe(field), 1, 0, "C" if i != 1 else "L", True)
            pdf.ln(6)
        pdf.ln(3)
        pdf.set_font("Helvetica", "I", 7)
        pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
        pdf.cell(w, 5, "All deterministic coordinate measurements within thresholds.", ln=True)
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        # ═══ REGULATOR FIX #3: LAW-100 inline footnote ═══
        pdf.ln(3)
        pdf.set_font("Helvetica", "I", 7)
        pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
        pdf.multi_cell(w, 4, "Note: LAW-100 (Bond Integrity RMSZ) is classified Deterministic for predicted structures "
                       "and Advisory for experimental coordinate regimes, reflecting resolution-dependent uncertainty. "
                       "See PIL-CAL-03 modality matrix for full classification rules.")
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])

        # ══════ PAGE 3: SECTION II - ADVISORY METRICS ══════
        pdf.add_page()

        # ═══ REGULATOR FIX #4: Distinct section header (different color) ═══
        pdf.section_divider("SECTION II: NON-SCORE-DRIVING CONTEXTUAL SIGNALS", 140, 100, 20)

        pdf.set_font("Helvetica", "I", 8)
        pdf.multi_cell(w, 4, "These metrics provide structural context but do NOT independently trigger VETO. "
                       "They are reported for transparency and expert review. Heuristic signals are statistical "
                       "proxies; Advisory signals are modality-reclassified deterministic laws.")
        pdf.ln(2)
        pdf.set_font("Helvetica", "B", 7)
        pdf.cell(w, 4, "Classification Definitions:", ln=True)
        pdf.set_font("Helvetica", "", 7)
        pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
        pdf.cell(w, 4, "DETERMINISTIC = Physically invariant coordinate measurement. Failure triggers VETO for predicted structures.", ln=True)
        pdf.cell(w, 4, "ADVISORY = Deterministic law reclassified to non-scoring for experimental data due to coordinate uncertainty (PIL-CAL-03).", ln=True)
        pdf.cell(w, 4, "HEURISTIC = Statistical proxy for structural plausibility. Does not trigger VETO. Contextual signal only.", ln=True)
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        pdf.ln(2)
        pdf.ln(2)

        pdf.ln(4)

        adv_laws = [l for l in laws if l.get('method') != 'deterministic']
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_fill_color(*PDF_COLORS["BG_ADV_TABLE"])
        for i, h in enumerate(headers):
            pdf.cell(widths[i], 8, h, 1, 0, "C", True)
        pdf.ln(8)
        pdf.set_font("Helvetica", "", 7)
        for l in adv_laws:
            cls = "Advisory" if "advisory" in l.get('method', '') else "Heuristic"
            display_status = l['status'] if l['status'] != 'VETO' else 'FAIL'
            for i, field in enumerate([l['law_id'], l['title'], cls,
                                       str(l['observed']), f"{l['operator']} {l['threshold']}", display_status]):
                pdf.cell(widths[i], 7, pdf.safe(field), 1, 0, "C" if i != 1 else "L")
            pdf.ln(7)

        # ══════ PAGE 4: CHARACTERIZATION ══════

        # ══════ STRUCTURAL VISUALIZATION PAGE ══════
        pdf.add_page()
        pdf.section_title("Structural Visualization")
        pdf.section_divider("3D CA-TRACE RENDER (CONFIDENCE-COLORED)", 20, 60, 120)
        pdf.ln(3)

        pdb_b64_data = payload.get("pdb_b64", "")
        ca_coords, ca_confidences = _extract_ca_atoms(pdb_b64_data)

        if len(ca_coords) > 500:
            step = max(1, len(ca_coords) // 500)
            ca_coords_plot = ca_coords[::step]
            ca_conf_plot = ca_confidences[::step]
        else:
            ca_coords_plot = ca_coords
            ca_conf_plot = ca_confidences

        if ca_coords_plot and len(ca_coords_plot) >= 3:
            try:
                _rtitle = "CA Trace - %s (%d residues)" % (prov.get("source", "Unknown"), len(ca_coords))
                render_png = generate_structure_render(ca_coords_plot, ca_conf_plot, _rtitle)
                import tempfile as _tf
                with _tf.NamedTemporaryFile(suffix=".png", delete=False) as _tmp:
                    _tmp.write(render_png)
                    _render_path = _tmp.name
                pdf.image(_render_path, x=15, w=180)
                import os as _os
                _os.unlink(_render_path)
                pdf.ln(5)
                pdf.set_font("Helvetica", "", 8)
                pdf.set_text_color(*PDF_COLORS["TEXT_MUTED"])
                _atom_msg = "Atoms plotted: %d CA | Colored by pLDDT/B-factor" % len(ca_coords)
                pdf.cell(w, 5, _atom_msg, ln=True, align="C")
                pdf.cell(w, 5, "Blue (>=90) | Light Blue (70-90) | Orange (50-70) | Red (<50)", ln=True, align="C")
                pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
            except Exception:
                pdf.set_font("Helvetica", "I", 9)
                pdf.cell(w, 8, "[Structure visualization could not be generated]", ln=True, align="C")
        else:
            pdf.set_font("Helvetica", "I", 9)
            pdf.cell(w, 8, "[Insufficient CA atoms for 3D visualization]", ln=True, align="C")

        # ── RESIDUE-LEVEL CONFIDENCE HEATMAP (below 3D render) ──
        if ca_confidences and len(ca_confidences) >= 3:
            try:
                import tempfile as _tf_hm
                import os as _os_hm
                _hm_title = "Per-Residue Confidence - %d residues" % len(ca_confidences)
                heatmap_png = generate_residue_heatmap_png(
                    ca_confidences, title=_hm_title,
                    figwidth=7.0, figheight=1.8, dpi=200,
                )
                with _tf_hm.NamedTemporaryFile(suffix=".png", delete=False) as _tmp_hm:
                    _tmp_hm.write(heatmap_png)
                    _hm_path = _tmp_hm.name
                if pdf.get_y() > 200:
                    pdf.add_page()
                    pdf.section_divider("RESIDUE-LEVEL CONFIDENCE DISTRIBUTION", 20, 60, 120)
                    pdf.ln(2)
                else:
                    pdf.ln(6)
                    pdf.set_font("Helvetica", "B", 9)
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.cell(w, 6, "RESIDUE-LEVEL CONFIDENCE DISTRIBUTION", ln=True)
                    pdf.ln(2)
                pdf.image(_hm_path, x=15, w=180)
                pdf.set_y(pdf.get_y() + 2)
                _os_hm.unlink(_hm_path)
                pdf.set_font("Helvetica", "", 7)
                pdf.set_text_color(*PDF_COLORS["TEXT_MUTED"])
                n_total = len(ca_confidences)
                n_low = sum(1 for s in ca_confidences if s < 50)
                n_med = sum(1 for s in ca_confidences if 50 <= s < 70)
                n_high = sum(1 for s in ca_confidences if 70 <= s < 90)
                n_vhigh = sum(1 for s in ca_confidences if s >= 90)
                pdf.cell(w, 4, "Very High %d (%.0f%%) | High %d (%.0f%%) | Medium %d (%.0f%%) | Low %d (%.0f%%)" % (
                    n_vhigh, n_vhigh/n_total*100, n_high, n_high/n_total*100,
                    n_med, n_med/n_total*100, n_low, n_low/n_total*100), ln=True, align="C")
                if n_low > 0:
                    low_regions = _find_low_regions(ca_confidences)
                    if low_regions:
                        pdf.set_text_color(*PDF_COLORS["VETO"])
                        pdf.set_font("Helvetica", "B", 7)
                        pdf.cell(w, 4, "Low-confidence regions: residues %s" % ", ".join(
                            ["%d-%d" % (s, e) for s, e in low_regions]), ln=True, align="C")
                pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                pdf.ln(2)
            except Exception:
                pdf.set_font("Helvetica", "I", 8)
                pdf.cell(w, 5, "[Residue heatmap could not be generated]", ln=True, align="C")
                pdf.ln(2)


        pdf.add_page()
        pdf.section_title("Structural Characterization")
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, f"Source Regime: {source_label}", ln=True)
        if ss.get("resolution"):
            pdf.cell(w, 8, f"Reported Resolution: {ss['resolution']} Angstroms", ln=True)
        pdf.ln(4)

        # ═══ REGULATOR FIX #5: Handle None/missing secondary structure ═══
        helix = ss.get('helix')
        sheet = ss.get('sheet')
        loop = ss.get('loop')

        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, "Secondary Structure Composition:", ln=True)
        pdf.set_font("Helvetica", "", 9)

        if helix is None and sheet is None and loop is None:
            if coverage < LAW_105_THRESHOLD:
                pdf.set_font("Helvetica", "I", 9)
                pdf.set_text_color(*PDF_COLORS["VETO"])
                pdf.multi_cell(w, 6, f"Secondary structure not computable: reliability coverage below {LAW_105_THRESHOLD}% "
                               f"(observed: {coverage}%). Insufficient high-confidence coordinates for "
                               "phi/psi-based assignment.")
                pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
            else:
                pdf.set_font("Helvetica", "I", 9)
                pdf.multi_cell(w, 6, "Secondary structure computation returned no data. "
                               "This may indicate a minimal peptide or non-standard chain topology.")
        else:
            h_val = f"{helix}%" if helix is not None else "N/A"
            s_val = f"{sheet}%" if sheet is not None else "N/A"
            l_val = f"{loop}%" if loop is not None else "N/A"
            pdf.cell(63, 10, f"Helix: {h_val}", 1, 0, "C")
            pdf.cell(63, 10, f"Sheet: {s_val}", 1, 0, "C")
            pdf.cell(64, 10, f"Loop: {l_val}", 1, 1, "C")

        pdf.ln(3)
        pdf.set_font("Helvetica", "I", 7)
        pdf.cell(w, 6, f"Scope: Computed over core residues only (n={ss.get('total_residues', 'N/A')}).", ln=True)

        # ══════ PAGE 5: PROVENANCE ══════

        # ══════ REFINEMENT RECOMMENDATIONS PAGE ══════
        if verdict_binary != "PASS" or coverage < LAW_105_THRESHOLD:
            pdf.add_page()
            pdf.section_title("Refinement Recommendations")
            pdf.section_divider("PHYSICS-FIRST REFINEMENT FOR DRUG DISCOVERY USE", 20, 60, 120)
            pdf.ln(3)

            pdf.set_font("Helvetica", "", 9)
            pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
            if coverage < LAW_105_THRESHOLD:
                pdf.multi_cell(w, 5,
                    "This structure has insufficient reliability coverage (%.1f%%) and is not yet "
                    "suitable for drug discovery without refinement. The following physics-first "
                    "methods are recommended to improve structural quality." % coverage)
            else:
                det_fail_list = [l for l in laws if l.get('method') == 'deterministic' and l['status'] != 'PASS']
                fail_names = ", ".join([l['law_id'] for l in det_fail_list])
                pdf.multi_cell(w, 5,
                    "This structure has deterministic violations (%s) that prevent certification. "
                    "The following refinement methods address the identified issues." % fail_names)
            pdf.ln(4)

            pdf.set_font("Helvetica", "B", 9)
            pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
            pdf.cell(w, 7, "Recommended Physics-First Refinement Methods (Prioritized)", ln=True)
            pdf.ln(2)

            ref_headers = ["Priority", "Method", "Expected Improvement", "Est. Time", "Recommendation"]
            ref_widths = [18, 40, 52, 22, 58]
            pdf.set_font("Helvetica", "B", 7)
            pdf.set_fill_color(*PDF_COLORS["SECTION_DET"])
            pdf.set_text_color(*PDF_COLORS["TEXT_WHITE"])
            for i, h in enumerate(ref_headers):
                pdf.cell(ref_widths[i], 7, h, 1, 0, "C", True)
            pdf.ln(7)

            ref_rows = [
                ["1", "Rosetta Fast Relax", "Fix rotamer outliers and clashes", "30-90 sec", "Strongly recommended"],
                ["2", "Short MD Equilibration", "Improve hydrophobic burial", "3-8 min", "Recommended"],
                ["3", "Targeted Loop Modeling", "Rebuild disordered regions", "2-5 min", "If loops are critical"],
                ["4", "Pocket Refinement", "Optimize binding site geometry", "8-15 min", "For drug design focus"],
            ]

            pdf.set_font("Helvetica", "", 7)
            pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
            for row_idx, row in enumerate(ref_rows):
                if row_idx % 2 == 0:
                    pdf.set_fill_color(*PDF_COLORS["BG_LIGHTER"])
                else:
                    pdf.set_fill_color(255, 255, 255)
                for i, field in enumerate(row):
                    pdf.cell(ref_widths[i], 6, pdf.safe(field), 1, 0,
                             "C" if i in [0, 3] else "L", True)
                pdf.ln(6)

            pdf.ln(5)

            # ── RESIDUE-LEVEL DIAGNOSTICS (LAW-125) ──────────────────────
            # Extract LAW-125 residue data if available
            law_125 = next((l for l in laws if l.get("law_id") == "LAW-125"), None)
            if law_125 and law_125.get("granularity") == "residue":
                rama_residues = law_125.get("failing_residues", [])
                if rama_residues:
                    pdf.set_font("Helvetica", "B", 9)
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.cell(w, 7, "Residue-Level Diagnostics", ln=True)
                    pdf.ln(1)
                    
                    pdf.set_font("Helvetica", "B", 8)
                    pdf.cell(w, 5, "LAW-125 Ramachandran Outliers:", ln=True)
                    pdf.set_font("Helvetica", "", 8)
                    
                    # Format: max 20 residues per line, institutional formatting
                    display_residues = rama_residues[:20]
                    residue_str = ", ".join(map(str, display_residues))
                    if len(rama_residues) > 20:
                        residue_str += " (+ %d more)" % (len(rama_residues) - 20)
                    
                    pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
                    pdf.multi_cell(w, 4,
                        "Detected %d Ramachandran outliers (%.2f%% of %d core residues): %s" % (
                            len(rama_residues),
                            law_125.get("observed", 0.0),
                            law_125.get("sample_size", law_125.get("sample", 0)),
                            residue_str
                        ))
                    pdf.ln(2)
                    
                    pdf.set_font("Helvetica", "I", 7)
                    pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
                    pdf.multi_cell(w, 4,
                        "Recommendation: Restrict Rosetta FastRelax to these specific residues "
                        "using coordinate constraints. Full-structure relaxation may introduce "
                        "new violations in currently compliant regions.")
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.ln(4)

            # ── RESIDUE-LEVEL DIAGNOSTICS (LAW-150) ──────────────────────
            law_150 = next((l for l in laws if l.get("law_id") == "LAW-150"), None)
            if law_150 and law_150.get("granularity") == "residue":
                rotamer_residues = law_150.get("failing_residues", [])
                if rotamer_residues:
                    pdf.set_font("Helvetica", "B", 8)
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.cell(w, 5, "LAW-150 Rotamer Outliers:", ln=True)
                    pdf.set_font("Helvetica", "", 8)

                    # Format: max 20 residues per line
                    display_residues = rotamer_residues[:20]
                    residue_str = ", ".join(map(str, display_residues))
                    if len(rotamer_residues) > 20:
                        residue_str += " (+ %d more)" % (len(rotamer_residues) - 20)

                    pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
                    pdf.multi_cell(w, 4,
                        "Detected %d rotamer outliers (%.2f%% of %d side chains): %s" % (
                            len(rotamer_residues),
                            law_150.get("observed", 0.0),
                            law_150.get("sample_size", law_150.get("sample", 0)),
                            residue_str
                        ))
                    pdf.ln(2)

                    pdf.set_font("Helvetica", "I", 7)
                    pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
                    pdf.multi_cell(w, 4,
                        "Recommendation: Apply Rosetta rotamer optimization with fa_dun scoreterm "
                        "upweighted to 2.0. These residues are outside Dunbrack rotamer library wells.")
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.ln(4)

            # ── RESIDUE-LEVEL DIAGNOSTICS (LAW-130) ──────────────────────
            law_130 = next((l for l in laws if l.get("law_id") == "LAW-130"), None)
            if law_130 and law_130.get("granularity") == "residue_pair":
                clash_pairs = law_130.get("failing_residue_pairs", [])
                if clash_pairs:
                    pdf.set_font("Helvetica", "B", 8)
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.cell(w, 5, "LAW-130 Steric Clashes:", ln=True)
                    pdf.set_font("Helvetica", "", 8)

                    # Format pairs: max 10 pairs displayed
                    display_pairs = clash_pairs[:10]
                    pairs_str = ", ".join([f"{p[0]}-{p[1]}" for p in display_pairs])
                    if len(clash_pairs) > 10:
                        pairs_str += " (+ %d more)" % (len(clash_pairs) - 10)

                    pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
                    pdf.multi_cell(w, 4,
                        "Detected %d clashing residue pairs (score: %.2f): %s" % (
                            len(clash_pairs),
                            law_130.get("observed", 0.0),
                            pairs_str
                        ))
                    pdf.ln(2)

                    pdf.set_font("Helvetica", "I", 7)
                    pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
                    pdf.multi_cell(w, 4,
                        "Recommendation: Apply Rosetta FastRelax with fa_rep upweighted to resolve steric overlaps.")
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.ln(4)

            # ── RESIDUE-LEVEL DIAGNOSTICS (LAW-135) ──────────────────────
            law_135 = next((l for l in laws if l.get("law_id") == "LAW-135"), None)
            if law_135 and law_135.get("granularity") == "residue":
                omega_residues = law_135.get("failing_residues", [])
                if omega_residues:
                    pdf.set_font("Helvetica", "B", 8)
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.cell(w, 5, "LAW-135 Omega Planarity Outliers:", ln=True)
                    pdf.set_font("Helvetica", "", 8)

                    display_residues = omega_residues[:20]
                    residue_str = ", ".join(map(str, display_residues))
                    if len(omega_residues) > 20:
                        residue_str += " (+ %d more)" % (len(omega_residues) - 20)

                    pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
                    pdf.multi_cell(w, 4,
                        "Detected %d omega outliers (%.2f%% of %d peptide bonds): %s" % (
                            len(omega_residues),
                            law_135.get("observed", 0.0),
                            law_135.get("sample_size", law_135.get("sample", 0)),
                            residue_str
                        ))
                    pdf.ln(2)

                    pdf.set_font("Helvetica", "I", 7)
                    pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
                    pdf.multi_cell(w, 4,
                        "Recommendation: Inspect cis-peptide bonds. Non-proline cis bonds are typically errors.")
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.ln(4)

            # ── RESIDUE-LEVEL DIAGNOSTICS (LAW-145) ──────────────────────
            law_145 = next((l for l in laws if l.get("law_id") == "LAW-145"), None)
            if law_145 and law_145.get("granularity") == "residue":
                chiral_residues = law_145.get("failing_residues", [])
                if chiral_residues:
                    pdf.set_font("Helvetica", "B", 8)
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.cell(w, 5, "LAW-145 Chirality Violations:", ln=True)
                    pdf.set_font("Helvetica", "", 8)

                    display_residues = chiral_residues[:20]
                    residue_str = ", ".join(map(str, display_residues))
                    if len(chiral_residues) > 20:
                        residue_str += " (+ %d more)" % (len(chiral_residues) - 20)

                    pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
                    pdf.multi_cell(w, 4,
                        "Detected %d D-amino acid chirality errors: %s" % (
                            len(chiral_residues),
                            residue_str
                        ))
                    pdf.ln(2)

                    pdf.set_font("Helvetica", "I", 7)
                    pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
                    pdf.multi_cell(w, 4,
                        "CRITICAL: Chirality errors require manual correction. These cannot be fixed by automated refinement.")
                    pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
                    pdf.ln(4)

            pdf.set_font("Helvetica", "B", 9)
            pdf.cell(w, 7, "Next Step", ln=True)
            pdf.set_font("Helvetica", "", 8)
            pdf.multi_cell(w, 5,
                "Run Rosetta Fast Relax on flagged regions, then re-upload the refined PDB "
                "to Toscanini for re-validation. After recommended refinement, this model is "
                "expected to reach sufficient reliability coverage and become suitable for "
                "virtual screening and lead optimization.")
            pdf.ln(3)

            pdf.set_font("Helvetica", "I", 7)
            pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
            pdf.multi_cell(w, 4,
                "Note: Refinement recommendations are generated based on detected physics "
                "violations and coverage analysis. Actual improvement depends on the specific "
                "structural context and refinement parameters used. All refined structures "
                "must be re-validated through the full Toscanini governance pipeline.")
            pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])


        pdf.add_page()
        pdf.section_title("Provenance & Authentication")
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, f"AUDIT_ID: {gov.get('audit_id', 'N/A')}", ln=True)
        pdf.cell(w, 8, f"TIMESTAMP_UTC: {gov.get('timestamp_utc', 'N/A')}", ln=True)
        pdf.cell(w, 8, f"CANON_HASH: {LAW_CANON_HASH}", ln=True)
        pdf.cell(w, 8, f"MATRIX_HASH: {compute_matrix_hash()}", ln=True)
        fp = gov.get("governance_fingerprint", {})
        if fp:
            pdf.cell(w, 8, f"MATRIX_SCHEMA: {fp.get('matrix_schema_version', 'N/A')}", ln=True)
            pdf.cell(w, 8, f"POLICY_REF: {fp.get('policy_ref', 'N/A')}", ln=True)
        pdf.cell(w, 8, f"ENGINE_VERSION: v{STATION_METADATA['version']}", ln=True)
        pdf.ln(5)

        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, "Coordinate SHA-256 Input Fingerprint:", ln=True)
        pdf.set_font("Courier", "", 8)
        pdf.multi_cell(w, 4, prov.get('hash', 'N/A'))

        pdf.ln(5)
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, "Strategic Math (PIL-MTH-12):", ln=True)
        pdf.set_font("Helvetica", "", 8)
        s6_val = math_data.get('s6', 'N/A')
        pdf.cell(w, 6, f"S6 (Deterministic Compliance): {s6_val}", ln=True)
        if s6_val == 0.0 or s6_val == 0:
            pdf.set_font("Helvetica", "I", 7)
            pdf.set_text_color(*PDF_COLORS["VETO"])
            pdf.multi_cell(w, 4, f"Note: S6 = 0.0 because coverage gate (LAW-105) failed. When reliability coverage < {LAW_105_THRESHOLD}%, the prioritization index is zeroed regardless of individual law compliance. This is by design: epistemic insufficiency prevents meaningful prioritization scoring.")
            pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
            pdf.set_font("Helvetica", "", 8)
        pdf.ln(3)
        pdf.set_x(15)
        pdf.set_font("Helvetica", "", 8)
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        pdf.cell(w, 6, f"W_arch (Architecture Weight): {math_data.get('w_arch', 'N/A')}", ln=True)
        pdf.set_x(15)
        pdf.cell(w, 6, f"M_S8 (NKG Penalty): {math_data.get('m_s8', 'N/A')}", ln=True)
        pdf.set_x(15)
        pdf.cell(w, 6, f"Formula: {BAYESIAN_FORMULA}", ln=True)

        return force_bytes(pdf.output(dest='S'))
    except Exception as e:
        import traceback
        return force_bytes(f"PDF GENERATION ERROR: {str(e)}\n{traceback.format_exc()}".encode())
