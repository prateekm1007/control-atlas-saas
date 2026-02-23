from fpdf import FPDF
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
        self.set_font("Helvetica", "B", 28)
        self.cell(badge_width, 22, verdict, border=0, ln=True, align="C", fill=True)

        # Score line below badge
        self.set_fill_color(*PDF_COLORS["BG_LIGHT"])
        self.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        self.set_font("Helvetica", "B", 11)
        self.cell(badge_width, 10,
                  f"Deterministic Score: {score}/100   |   Coverage: {coverage}%",
                  border=0, ln=True, align="C", fill=True)
        self.ln(4)

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

        # VERDICT BADGE — Large centered institutional panel
        det_score = v.get("deterministic_score", 0)
        pdf.verdict_badge(verdict_binary, det_score, coverage)

        # METRIC CARDS — Three columns, no borders
        card_w = 63

        # Labels row
        pdf.set_fill_color(*PDF_COLORS["BG_LIGHT"])
        pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
        pdf.set_font("Helvetica", "", 8)
        pdf.cell(card_w, 6, "DETERMINISTIC COMPLIANCE", border=0, align="C", fill=True)
        pdf.cell(card_w, 6, "RELIABILITY COVERAGE", border=0, align="C", fill=True)
        pdf.cell(card_w, 6, "PRIORITIZATION INDEX", border=0, ln=True, align="C", fill=True)

        # Values row
        pdf.set_fill_color(255, 255, 255)
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
        pdf.set_font("Helvetica", "B", 14)
        pdf.cell(card_w, 10, f"{v.get('det_passed', 0)}/{v.get('det_total', 12)}", border=1, align="C", fill=True)
        pdf.cell(card_w, 10, f"{coverage}%", border=1, align="C", fill=True)
        pdf.cell(card_w, 10, f"{payload.get('tier3', {}).get('probability', 0)}%", border=1, ln=True, align="C", fill=True)

        pdf.ln(3)

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

        # ══════ PAGE 2: SECTION I — DETERMINISTIC GOVERNANCE ══════
        pdf.add_page()

        # ═══ REGULATOR FIX #4: Distinct section header ═══
        pdf.section_divider("SECTION I: DETERMINISTIC GOVERNANCE (SCORE-DRIVING)", 20, 60, 120)

        pdf.set_font("Helvetica", "I", 8)
        pdf.multi_cell(w, 4, "Physical Refusal Gate. Failure of any deterministic law triggers VETO for predicted models. "
                       "For experimental structures, modality-specific reclassifications apply per PIL-CAL-03.")
        pdf.ln(2)
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_text_color(*PDF_COLORS["ALERT"])
        # CANONICAL ADJUDICATION HIERARCHY — DO NOT DUPLICATE ANYWHERE
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
        pdf.cell(w, 8, "DETERMINISTIC COMPLIANCE LEDGER", ln=True)
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_fill_color(*PDF_COLORS["BG_PASS_TABLE"])
        for i, h in enumerate(headers):
            pdf.cell(widths[i], 8, h, 1, 0, "C", True)
        pdf.ln(8)
        pdf.set_font("Helvetica", "", 7)
        for l in passes:
            for i, field in enumerate([l['law_id'], l['title'], "Deterministic",
                                       str(l['observed']), f"{l['operator']} {l['threshold']}", l['status']]):
                pdf.cell(widths[i], 7, pdf.safe(field), 1, 0, "C" if i != 1 else "L")
            pdf.ln(7)

        # ═══ REGULATOR FIX #3: LAW-100 inline footnote ═══
        pdf.ln(3)
        pdf.set_font("Helvetica", "I", 7)
        pdf.set_text_color(*PDF_COLORS["TEXT_SUBTLE"])
        pdf.multi_cell(w, 4, "Note: LAW-100 (Bond Integrity RMSZ) is classified Deterministic for predicted structures "
                       "and Advisory for experimental coordinate regimes, reflecting resolution-dependent uncertainty. "
                       "See PIL-CAL-03 modality matrix for full classification rules.")
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])

        # ══════ PAGE 3: SECTION II — ADVISORY METRICS ══════
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
        pdf.set_font("Helvetica", "B", 7)
        pdf.cell(w, 4, "Classification Definitions:", ln=True)
        pdf.set_font("Helvetica", "", 7)
        pdf.set_text_color(*PDF_COLORS["TEXT_SECONDARY"])
        pdf.cell(w, 4, "DETERMINISTIC = Physically invariant coordinate measurement. Failure triggers VETO for predicted structures.", ln=True)
        pdf.cell(w, 4, "ADVISORY = Deterministic law reclassified to non-scoring for experimental data due to coordinate uncertainty (PIL-CAL-03).", ln=True)
        pdf.cell(w, 4, "HEURISTIC = Statistical proxy for structural plausibility. Does not trigger VETO. Contextual signal only.", ln=True)
        pdf.set_text_color(*PDF_COLORS["TEXT_PRIMARY"])
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
            for i, field in enumerate([l['law_id'], l['title'], cls,
                                       str(l['observed']), f"{l['operator']} {l['threshold']}", l['status']]):
                pdf.cell(widths[i], 7, pdf.safe(field), 1, 0, "C" if i != 1 else "L")
            pdf.ln(7)

        # ══════ PAGE 4: CHARACTERIZATION ══════
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
        pdf.cell(w, 6, f"W_arch (Architecture Weight): {math_data.get('w_arch', 'N/A')}", ln=True)
        pdf.cell(w, 6, f"M_S8 (NKG Penalty): {math_data.get('m_s8', 'N/A')}", ln=True)
        pdf.cell(w, 6, f"Formula: {BAYESIAN_FORMULA}", ln=True)

        return force_bytes(pdf.output(dest='S'))
    except Exception as e:
        import traceback
        return force_bytes(f"PDF GENERATION ERROR: {str(e)}\n{traceback.format_exc()}".encode())
