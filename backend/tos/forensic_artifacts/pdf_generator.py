from fpdf import FPDF
from datetime import datetime
from ..utils.type_guards import force_bytes
from ..governance.station_sop import (
    STATION_METADATA, BAYESIAN_FORMULA, LAW_CANON_HASH
)

class ToscaniniDossier(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 8)
        self.set_text_color(80, 80, 80)
        self.cell(w=0, h=10, txt=f"TOSCANINI STRUCTURAL GOVERNANCE DOSSIER v{STATION_METADATA['version']}", border=0, align="L")
        self.ln(10)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 7)
        self.set_text_color(128, 128, 128)
        self.cell(0, 10, f"Engine: v{STATION_METADATA['version']} // SHA-256 Validated // Page {self.page_no()}", align="R")

    def section_title(self, text):
        self.set_font("Helvetica", "B", 11)
        self.set_text_color(0, 0, 0)
        self.cell(0, 8, text.upper(), ln=True)
        self.set_draw_color(0, 0, 0)
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(4)

    def safe(self, text):
        return str(text).encode('latin-1', 'replace').decode('latin-1')

def generate_v21_dossier(payload):
    try:
        pdf = ToscaniniDossier()
        pdf.set_auto_page_break(True, margin=15)
        w = 190
        v = payload.get('verdict', {})
        ss = payload.get("characterization", {})
        prov = payload.get("provenance", {})
        conf_meta = payload.get("confidence_meta", {})
        gov = payload.get("governance", {})
        math_data = payload.get("strategic_math", {})
        laws = payload.get('tier1', {}).get('laws', [])
        
        is_experimental = ss.get("source_type") != "predicted"

        # ══════ PAGE 1: EXECUTIVE DETERMINATION ══════
        pdf.add_page()
        pdf.section_title("Executive Determination")
        
        # Source Attribution Header
        source_label = f"EXPERIMENTAL ({ss.get('source_type', 'N/A').upper()})" if is_experimental else "ML-PREDICTED (ALPHAFOLD)"
        pdf.set_font("Helvetica", "B", 10)
        pdf.set_fill_color(240, 240, 240)
        pdf.cell(w, 10, f"DATA REGIME: {source_label}", border=1, ln=True, fill=True, align="C")
        pdf.ln(5)

        decision_data = [
            ["Final Verdict", v.get("binary", "ERROR")],
            ["Deterministic Compliance", f"{v.get('det_passed', 0)} / {v.get('det_total', 11)}"],
            ["Reliability Coverage", f"{v.get('coverage_pct', 0)} %"],
            ["Prioritization Index", f"{payload.get('tier3', {}).get('probability', 0)} %"]
        ]

        for row in decision_data:
            pdf.set_font("Helvetica", "B", 9)
            pdf.set_fill_color(250, 250, 250)
            pdf.cell(80, 8, row[0], border=1, fill=True)
            pdf.set_font("Helvetica", "", 9)
            if row[1] == "VETO": pdf.set_text_color(170, 35, 35)
            elif row[1] == "PASS": pdf.set_text_color(15, 110, 55)
            pdf.cell(110, 8, row[1], border=1, ln=True)
            pdf.set_text_color(0, 0, 0)

        pdf.ln(5)
        # 🛡️ Regulator-Grade Interpretation Clause
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "INTERPRETATION OF VERDICT CLASSIFICATION", ln=True)
        pdf.set_font("Helvetica", "", 8)
        
        pdf.set_font("Helvetica", "B", 8); pdf.write(5, "1. PASS: "); pdf.set_font("Helvetica", "", 8)
        pdf.write(5, "Physically compliant with defined thresholds. Coverage >= 70%.\n")
        
        pdf.set_font("Helvetica", "B", 8); pdf.write(5, "2. VETO: "); pdf.set_font("Helvetica", "", 8)
        pdf.write(5, "Physically incompatible or internally inconsistent condition detected. Violation of non-negotiable physical canon.\n")
        
        pdf.set_font("Helvetica", "B", 8); pdf.write(5, "3. INDETERMINATE: "); pdf.set_font("Helvetica", "", 8)
        pdf.write(5, "Engine asserts epistemic insufficiency to certify compliance. Reflects coverage < 70% (LAW-105) or low-confidence data. Not a declaration of failure.\n")
        
        pdf.ln(4)
        pdf.set_font("Helvetica", "I", 8)
        pdf.multi_cell(w, 4, "Note: For experimental structures, geometric metrics (LAW-100) are classified as Advisory to account for resolution-based coordinate uncertainty.")

        # ══════ PAGE 2: SECTION I — DETERMINISTIC GOVERNANCE ══════
        pdf.add_page()
        pdf.section_title("Section I: Deterministic Governance (Score-Driving)")
        pdf.set_font("Helvetica", "I", 8)
        pdf.multi_cell(w, 4, "Physical 'Refusal Gate'. Failure of any law in this section triggers VETO for predicted models. Coverage < 70% invalidates adjudication per LAW-105.")
        pdf.ln(4)

        # Split: Violations vs Passes (Removing Redundancy)
        det_laws = [l for l in laws if l.get('method') == 'deterministic']
        violations = [l for l in det_laws if l['status'] != 'PASS']
        passes = [l for l in det_laws if l['status'] == 'PASS']

        headers = ["ID", "Metric", "Class", "Observed", "Threshold", "Status"]
        widths = [15, 50, 30, 35, 35, 25]

        if violations:
            pdf.set_font("Helvetica", "B", 9); pdf.set_text_color(170, 35, 35)
            pdf.cell(w, 8, f"DETECTION: {len(violations)} DETERMINISTIC VIOLATIONS", ln=True)
            pdf.set_text_color(0,0,0); pdf.set_font("Helvetica", "B", 8); pdf.set_fill_color(245, 230, 230)
            for i, h in enumerate(headers): pdf.cell(widths[i], 8, h, 1, 0, "C", True)
            pdf.ln(8); pdf.set_font("Helvetica", "", 7)
            for l in violations:
                for i, field in enumerate([l['law_id'], l['title'], "Deterministic", str(l['observed']), f"{l['operator']} {l['threshold']}", l['status']]):
                    pdf.cell(widths[i], 7, pdf.safe(field), 1, 0, "C" if i != 1 else "L")
                pdf.ln(7)
            pdf.ln(5)

        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, "DETERMINISTIC COMPLIANCE LEDGER", ln=True)
        pdf.set_font("Helvetica", "B", 8); pdf.set_fill_color(230, 245, 230)
        for i, h in enumerate(headers): pdf.cell(widths[i], 8, h, 1, 0, "C", True)
        pdf.ln(8); pdf.set_font("Helvetica", "", 7)
        for l in passes:
            for i, field in enumerate([l['law_id'], l['title'], "Deterministic", str(l['observed']), f"{l['operator']} {l['threshold']}", l['status']]):
                pdf.cell(widths[i], 7, pdf.safe(field), 1, 0, "C" if i != 1 else "L")
            pdf.ln(7)

        # ══════ PAGE 3: SECTION II — ADVISORY METRICS ══════
        pdf.add_page()
        pdf.section_title("Section II: Advisory & Clinical Signals")
        pdf.set_font("Helvetica", "I", 8)
        pdf.multi_cell(w, 4, "Contextual signals. These metrics do not independently trigger VETO. LAW-100 (RMSZ) is advisory for experimental data regimes.")
        pdf.ln(4)

        adv_laws = [l for l in laws if l.get('method') != 'deterministic']
        pdf.set_font("Helvetica", "B", 8); pdf.set_fill_color(255, 250, 230)
        for i, h in enumerate(headers): pdf.cell(widths[i], 8, h, 1, 0, "C", True)
        pdf.ln(8); pdf.set_font("Helvetica", "", 7)
        for l in adv_laws:
            cls = "Advisory" if "advisory" in l['method'] else "Heuristic"
            for i, field in enumerate([l['law_id'], l['title'], cls, str(l['observed']), f"{l['operator']} {l['threshold']}", l['status']]):
                pdf.cell(widths[i], 7, pdf.safe(field), 1, 0, "C" if i != 1 else "L")
            pdf.ln(7)

        # ══════ PAGE 4: CHARACTERIZATION ══════
        pdf.add_page()
        pdf.section_title("Structural Characterization")
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, f"Source Regime: {source_label}", ln=True)
        if ss.get("resolution"): pdf.cell(w, 8, f"Reported Resolution: {ss['resolution']} Angstroms", ln=True)
        pdf.ln(4)

        pdf.cell(w, 8, "Secondary Structure Composition:", ln=True)
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(63, 10, f"Helix: {ss.get('helix')}%", 1, 0, "C")
        pdf.cell(63, 10, f"Sheet: {ss.get('sheet')}%", 1, 0, "C")
        pdf.cell(64, 10, f"Loop: {ss.get('loop')}%", 1, 1, "C")
        pdf.set_font("Helvetica", "I", 7)
        pdf.cell(w, 6, f"Scope Note: Computed over core residues only (n={ss.get('total_residues')}).", ln=True)

        # ══════ PAGE 5: PROVENANCE ══════
        pdf.add_page()
        pdf.section_title("Provenance & Authentication")
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, f"AUDIT_ID: {gov.get('audit_id')}", ln=True)
        pdf.cell(w, 8, f"TIMESTAMP_UTC: {gov.get('timestamp_utc')}", ln=True)
        pdf.cell(w, 8, f"CANON_HASH: {LAW_CANON_HASH}", ln=True)
        pdf.ln(5)
        pdf.cell(w, 8, "Coordinate SHA-256 Input Fingerprint:", ln=True)
        pdf.set_font("Courier", "", 8); pdf.multi_cell(w, 4, prov.get('hash', 'N/A'))

        return force_bytes(pdf.output(dest='S'))
    except Exception as e:
        import traceback
        return force_bytes(f"PDF GENERATION ERROR: {str(e)}\n{traceback.format_exc()}".encode())
