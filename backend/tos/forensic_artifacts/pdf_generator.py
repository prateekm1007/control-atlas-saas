from fpdf import FPDF
from datetime import datetime
from ..utils.type_guards import sanitize_notary_text, force_bytes
from ..governance.station_sop import (
    LAW_CANON, SCORE_DEFINITIONS, STATION_METADATA, BAYESIAN_FORMULA,
    LAW_CATEGORIES, CATEGORY_DISPLAY, get_laws_by_category, LAW_CANON_HASH
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
        # 🛡️ PIL-NOT-18: Institutional Footer
        self.cell(0, 10, f"Engine: v{STATION_METADATA['version']} // SHA-256 Validated // Page {self.page_no()}", align="R")
        
    def section_title(self, text):
        self.set_font("Helvetica", "B", 11)
        self.set_text_color(0, 0, 0)
        self.cell(0, 8, text.upper(), ln=True)
        self.set_draw_color(0, 0, 0)
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(4)

    def safe(self, text):
        # 🛡️ PIL-HYG-21: Support for IUPAC symbols via Latin-1 mapping
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
        math = payload.get("strategic_math", {})
        laws = payload.get('tier1', {}).get('laws', [])
        
        # ══════ PAGE 1: EXECUTIVE DETERMINATION ══════
        pdf.add_page()
        pdf.section_title("Executive Determination")
        
        pdf.set_font("Helvetica", "B", 9)
        decision_data = [
            ["Final Verdict", v.get("binary", "ERROR")],
            ["Deterministic Failures", f"{v.get('det_total', 12) - v.get('det_passed', 0)} / {v.get('det_total', 12)}"],
            ["Heuristic Flags", f"{v.get('heur_total', 3) - v.get('heur_passed', 0)} / {v.get('heur_total', 3)}"],
            ["Reliability Coverage", f"{v.get('coverage_pct', 0)} %"],
            ["Mean pLDDT", f"{round(v.get('confidence_score', 0), 2)}"],
            ["Composite Prioritization Index", f"{payload.get('tier3', {}).get('probability', 0)} %"]
        ]
        
        for row in decision_data:
            pdf.set_font("Helvetica", "B", 9)
            pdf.set_fill_color(245, 245, 245)
            pdf.cell(80, 8, row[0], border=1, fill=True)
            pdf.set_font("Helvetica", "", 9)
            # High-visibility verdict highlight
            if row[0] == "Final Verdict" and row[1] == "VETO": pdf.set_text_color(170, 35, 35)
            elif row[0] == "Final Verdict" and row[1] == "PASS": pdf.set_text_color(15, 110, 55)
            pdf.cell(110, 8, row[1], border=1, ln=True)
            pdf.set_text_color(0, 0, 0)
        
        pdf.ln(8)
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 7, "Governance Recommendation:", ln=True)
        pdf.set_font("Helvetica", "I", 9)
        rec = "Insufficient coverage gate. Re-evaluate with experimental data."
        if v.get("binary") == "PASS": rec = "Structure compliant with physical invariants. Eligible for lead selection."
        elif v.get("binary") == "VETO": rec = "Deterministic physics violated. Do not allocate downstream resources."
        pdf.multi_cell(w, 5, rec)

        # ══════ PAGE 2: DETERMINISTIC FAILURES ══════
        pdf.add_page()
        pdf.section_title("Critical Invariant Violations (VETO Gate)")
        failures = [l for l in laws if l['status'] in ("FAIL", "VETO") and l['method'] == "deterministic" and l['law_id'] != "LAW-105"]
        
        if not failures:
            pdf.set_font("Helvetica", "B", 10)
            pdf.set_text_color(0, 100, 0)
            pdf.cell(w, 15, "NO DETERMINISTIC INVARIANT VIOLATIONS DETECTED", align="C", ln=True)
        else:
            pdf.set_font("Helvetica", "B", 8)
            pdf.set_fill_color(240, 230, 230)
            headers = ["ID", "Observed", "Threshold", "Excess", "Units"]
            widths = [20, 45, 45, 40, 40]
            for i, h in enumerate(headers): pdf.cell(widths[i], 8, h, 1, 0, "C", True)
            pdf.ln(8)
            pdf.set_font("Helvetica", "", 8)
            pdf.set_text_color(170, 35, 35)
            for l in failures:
                pdf.cell(20, 8, l['law_id'], 1, 0, "C")
                pdf.cell(45, 8, str(l['observed']), 1, 0, "C")
                pdf.cell(45, 8, f"{l['operator']} {l['threshold']}", 1, 0, "C")
                pdf.cell(40, 8, l['deviation'], 1, 0, "C")
                pdf.cell(40, 8, pdf.safe(l['units']), 1, 1, "C")
        pdf.set_text_color(0, 0, 0)

        # ══════ PAGE 3: DETERMINISTIC PASS TABLE ══════
        pdf.add_page()
        pdf.section_title("Deterministic Invariants (Compliance Table)")
        passes = [l for l in laws if l['status'] == "PASS" and l['method'] == "deterministic"]
        
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_fill_color(235, 245, 235)
        headers = ["ID", "Metric", "Observed", "Threshold", "Dev", "Units"]
        widths = [18, 52, 35, 35, 25, 25]
        for i, h in enumerate(headers): pdf.cell(widths[i], 8, h, 1, 0, "C", True)
        pdf.ln(8)
        
        pdf.set_font("Helvetica", "", 7)
        for l in passes:
            pdf.cell(18, 7, l['law_id'], 1, 0, "C")
            pdf.cell(52, 7, l['title'], 1)
            pdf.cell(35, 7, str(l['observed']), 1, 0, "C")
            pdf.cell(35, 7, f"{l['operator']} {l['threshold']}", 1, 0, "C")
            pdf.cell(25, 7, l['deviation'], 1, 0, "C")
            pdf.cell(25, 7, pdf.safe(l['units']), 1, 1, "C")
            
        pdf.ln(5)
        pdf.set_font("Helvetica", "I", 7)
        pdf.multi_cell(w, 4, "Reference Standard: Engh-Huber 1991. ML-Predicted structures evaluated against same invariant reference standards.")

        # ══════ PAGE 4: HEURISTIC SIGNALS ══════
        pdf.add_page()
        pdf.section_title("Heuristic Advisory Signals")
        heuristics = [l for l in laws if l['method'] == "heuristic"]
        
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_fill_color(245, 245, 240)
        for i, h in enumerate(headers): pdf.cell(widths[i], 8, h, 1, 0, "C", True)
        pdf.ln(8)
        
        pdf.set_font("Helvetica", "", 7)
        for l in heuristics:
            pdf.cell(18, 7, l['law_id'], 1, 0, "C")
            pdf.cell(52, 7, l['title'], 1)
            pdf.cell(35, 7, str(l['observed']), 1, 0, "C")
            pdf.cell(35, 7, f"{l['operator']} {l['threshold']}", 1, 0, "C")
            pdf.cell(25, 7, l['deviation'], 1, 0, "C")
            pdf.cell(25, 7, pdf.safe(l['units']), 1, 1, "C")

        # ══════ PAGE 5: STRUCTURAL CHARACTERIZATION ══════
        pdf.add_page()
        pdf.section_title("Structural Characterization")
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, f"Secondary Structure (Method: {ss.get('method')})", ln=True)
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(63, 10, f"Alpha Helix: {ss.get('helix')}%", 1, 0, "C")
        pdf.cell(63, 10, f"Beta Sheet: {ss.get('sheet')}%", 1, 0, "C")
        pdf.cell(64, 10, f"Loop/Coil: {ss.get('loop')}%", 1, 1, "C")
        pdf.set_font("Helvetica", "I", 7)
        pdf.multi_cell(w, 5, f"Analysis Denominator: Percentages computed over high-confidence core residues only (n={ss.get('core_evaluated')}).")
        pdf.ln(5)
        
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "Confidence Distribution (pLDDT)", ln=True)
        dist = conf_meta.get("distribution", {})
        if dist:
            pdf.cell(63, 10, f"Core (>=70): {dist.get('high_conf_pct')}%", 1, 0, "C")
            pdf.cell(63, 10, f"Boundary (50-70): {dist.get('med_conf_pct')}%", 1, 0, "C")
            pdf.cell(64, 10, f"Fringe (<50): {dist.get('low_conf_pct')}%", 1, 1, "C")

        # ══════ PAGE 6: PRIORITIZATION MODEL ══════
        pdf.add_page()
        pdf.section_title("Prioritization Model Expansion")
        pdf.set_font("Courier", "B", 10)
        pdf.cell(w, 10, f"FORMULA: {BAYESIAN_FORMULA}", ln=True, align="C")
        pdf.ln(5)
        
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "Numeric Expansion:", ln=True)
        pdf.set_font("Courier", "", 10)
        pdf.cell(w, 8, f"S6 (Integrity Factor) = {math.get('s6', 0.0)}", ln=True)
        pdf.cell(w, 8, f"W_arch ({math.get('architecture')}) = {math.get('w_arch', 1.0)}", ln=True)
        pdf.cell(w, 8, f"M_S8 (Negative Penalty) = {math.get('m_s8', 0.0)}", ln=True)
        pdf.ln(5)
        
        if v.get("suppression_reason"):
            pdf.set_font("Helvetica", "B", 10)
            pdf.set_text_color(170, 35, 35)
            pdf.cell(w, 10, f"OVERRIDE: {v['suppression_reason']}", border=1, ln=True, align="C")
            pdf.set_text_color(0, 0, 0)

        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(w, 15, f"COMPOSITE PRIORITIZATION SCORE (Non-Probabilistic): {payload.get('tier3', {}).get('probability', 0)}%", border=1, align="C")

        # ══════ PAGE 7: PROVENANCE & AUTHENTICATION ══════
        pdf.add_page()
        pdf.section_title("Provenance & Authentication")
        pdf.set_font("Courier", "B", 9)
        pdf.cell(w, 8, f"AUDIT_ID:       {gov.get('audit_id')}", ln=True)
        pdf.cell(w, 8, f"SOURCE:         {prov.get('source')}", ln=True)
        pdf.cell(w, 8, f"TIMESTAMP_UTC:  {gov.get('timestamp_utc')}", ln=True)
        pdf.cell(w, 8, f"TOTAL_RESIDUES: {ss.get('total_residues')}", ln=True)
        pdf.cell(w, 8, f"TOTAL_ATOMS:    {ss.get('total_atoms')}", ln=True)
        pdf.ln(5)
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "Input Authentication (SHA-256):", ln=True)
        pdf.set_font("Courier", "", 8)
        pdf.multi_cell(w, 4, prov.get('hash'))
        pdf.ln(5)
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "Canon Fingerprint (Truncated SHA-256, 16 bytes):", ln=True)
        pdf.set_font("Courier", "", 8)
        pdf.cell(w, 8, LAW_CANON_HASH, ln=True)

        return force_bytes(pdf.output(dest='S'))
    except Exception as e:
        err = FPDF()
        err.add_page()
        err.set_font("Helvetica", "B", 12)
        err.cell(0, 10, f"CRITICAL REPORTING ERROR: {str(e)}", ln=True)
        return force_bytes(err.output(dest='S'))

