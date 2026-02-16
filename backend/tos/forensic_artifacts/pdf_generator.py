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

        # Classify laws
        det_laws = [l for l in laws if l.get('method') == 'deterministic']
        adv_exp_laws = [l for l in laws if l.get('method') == 'advisory_experimental']
        heur_laws = [l for l in laws if l.get('method') == 'heuristic']

        det_passed = sum(1 for l in det_laws if l['status'] == 'PASS')
        det_total = len(det_laws)
        has_advisory = len(adv_exp_laws) > 0
        source_type = ss.get("source_type", "predicted")

        # ══════ PAGE 1: EXECUTIVE DETERMINATION ══════
        pdf.add_page()
        pdf.section_title("Executive Determination")

        # Source badge
        if source_type == "experimental":
            pdf.set_font("Helvetica", "B", 9)
            pdf.set_fill_color(230, 240, 255)
            pdf.cell(w, 8, f"SOURCE: EXPERIMENTAL ({ss.get('resolution', 'N/A')} A resolution)", border=1, ln=True, fill=True, align="C")
            pdf.ln(2)
        else:
            pdf.set_font("Helvetica", "B", 9)
            pdf.set_fill_color(255, 245, 230)
            pdf.cell(w, 8, "SOURCE: PREDICTED (AlphaFold)", border=1, ln=True, fill=True, align="C")
            pdf.ln(2)

        det_label = f"{det_passed} / {det_total}"
        if has_advisory:
            det_label += " + 1 Advisory"

        decision_data = [
            ["Final Verdict", v.get("binary", "ERROR")],
            ["Deterministic Compliance", det_label],
            ["Heuristic Flags", f"{v.get('heur_total', 3) - v.get('heur_passed', 0)} / {v.get('heur_total', 3)}"],
            ["Reliability Coverage", f"{v.get('coverage_pct', 0)} %"],
        ]

        if source_type == "experimental":
            decision_data.append(["Resolution", f"{ss.get('resolution', 'N/A')} A"])
            decision_data.append(["Mean B-factor", f"{conf_meta.get('mean', 'N/A')}"])
        else:
            decision_data.append(["Mean pLDDT", f"{round(v.get('confidence_score', 0), 2)}"])

        decision_data.append(["Prioritization Index", f"{payload.get('tier3', {}).get('probability', 0)} %"])

        for row in decision_data:
            pdf.set_font("Helvetica", "B", 9)
            pdf.set_fill_color(245, 245, 245)
            pdf.cell(80, 8, row[0], border=1, fill=True)
            pdf.set_font("Helvetica", "", 9)
            if row[0] == "Final Verdict" and row[1] == "VETO":
                pdf.set_text_color(170, 35, 35)
            elif row[0] == "Final Verdict" and row[1] == "PASS":
                pdf.set_text_color(15, 110, 55)
            pdf.cell(110, 8, row[1], border=1, ln=True)
            pdf.set_text_color(0, 0, 0)

        pdf.ln(5)
        pdf.set_font("Helvetica", "I", 9)
        rec = "Insufficient coverage. Re-evaluate with experimental data."
        if v.get("binary") == "PASS":
            rec = "Structure compliant with physical invariants. Eligible for lead selection."
        elif v.get("binary") == "VETO":
            rec = "Deterministic physics violated. Do not allocate downstream resources."
        pdf.multi_cell(w, 5, rec)

        # ══════ PAGE 2: SECTION I — DETERMINISTIC GOVERNANCE ══════
        pdf.add_page()
        pdf.section_title("Section I: Deterministic Governance (Score-Driving)")

        pdf.set_font("Helvetica", "I", 8)
        pdf.multi_cell(w, 4, f"These {det_total} laws are physically deterministic. Each must PASS for the structure to receive a PASS verdict. Failure of any single law triggers VETO.")
        pdf.ln(3)

        # Failures first
        det_failures = [l for l in det_laws if l['status'] in ("FAIL", "VETO")]
        if det_failures:
            pdf.set_font("Helvetica", "B", 9)
            pdf.set_text_color(170, 35, 35)
            pdf.cell(w, 8, f"DETERMINISTIC VIOLATIONS: {len(det_failures)}", ln=True)
            pdf.set_text_color(0, 0, 0)

            pdf.set_font("Helvetica", "B", 8)
            pdf.set_fill_color(240, 230, 230)
            hdrs = ["ID", "Metric", "Observed", "Threshold", "Dev", "Units"]
            wds = [18, 52, 35, 35, 25, 25]
            for i, h in enumerate(hdrs):
                pdf.cell(wds[i], 8, h, 1, 0, "C", True)
            pdf.ln(8)
            pdf.set_font("Helvetica", "", 7)
            pdf.set_text_color(170, 35, 35)
            for l in det_failures:
                pdf.cell(18, 7, l['law_id'], 1, 0, "C")
                pdf.cell(52, 7, l['title'], 1)
                pdf.cell(35, 7, str(l['observed']), 1, 0, "C")
                pdf.cell(35, 7, f"{l['operator']} {l['threshold']}", 1, 0, "C")
                pdf.cell(25, 7, l['deviation'], 1, 0, "C")
                pdf.cell(25, 7, pdf.safe(l['units']), 1, 1, "C")
            pdf.set_text_color(0, 0, 0)
            pdf.ln(5)

        # All deterministic (compliance table)
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, f"Full Deterministic Compliance ({det_passed}/{det_total} PASS)", ln=True)

        pdf.set_font("Helvetica", "B", 8)
        pdf.set_fill_color(235, 245, 235)
        hdrs = ["ID", "Metric", "Observed", "Threshold", "Dev", "Units"]
        wds = [18, 52, 35, 35, 25, 25]
        for i, h in enumerate(hdrs):
            pdf.cell(wds[i], 8, h, 1, 0, "C", True)
        pdf.ln(8)
        pdf.set_font("Helvetica", "", 7)
        for l in det_laws:
            if l['status'] == 'PASS':
                pdf.set_text_color(0, 80, 0)
            else:
                pdf.set_text_color(170, 35, 35)
            pdf.cell(18, 7, l['law_id'], 1, 0, "C")
            pdf.cell(52, 7, l['title'], 1)
            pdf.cell(35, 7, str(l['observed']), 1, 0, "C")
            pdf.cell(35, 7, f"{l['operator']} {l['threshold']}", 1, 0, "C")
            pdf.cell(25, 7, l['deviation'], 1, 0, "C")
            pdf.cell(25, 7, pdf.safe(l['units']), 1, 1, "C")
        pdf.set_text_color(0, 0, 0)

        pdf.ln(3)
        pdf.set_font("Helvetica", "I", 7)
        pdf.multi_cell(w, 4, "Reference: Engh & Huber (1991), Acta Cryst. A47. ML-predicted structures evaluated against identical invariant thresholds.")

        # ══════ PAGE 3: SECTION II — ADVISORY METRICS ══════
        pdf.add_page()
        pdf.section_title("Section II: Advisory Metrics")

        if adv_exp_laws:
            pdf.set_font("Helvetica", "B", 9)
            pdf.set_fill_color(255, 248, 220)
            pdf.cell(w, 8, "ADVISORY (EXPERIMENTAL) — Not scored, reported for transparency", border=1, ln=True, fill=True, align="C")
            pdf.ln(3)

            pdf.set_font("Helvetica", "I", 8)
            pdf.multi_cell(w, 4,
                "Advisory laws are excluded from deterministic adjudication for experimental structures. "
                "Engh-Huber bond length sigma values describe restraint targets in small-molecule crystals, "
                "not coordinate uncertainty at macromolecular resolution (1.5-3.0 A). "
                "RMSZ is reported for comparative reference only."
            )
            pdf.ln(3)

            pdf.set_font("Helvetica", "B", 8)
            pdf.set_fill_color(255, 248, 220)
            hdrs_adv = ["ID", "Metric", "Observed (RMSZ)", "Classification"]
            wds_adv = [25, 60, 55, 50]
            for i, h in enumerate(hdrs_adv):
                pdf.cell(wds_adv[i], 8, h, 1, 0, "C", True)
            pdf.ln(8)
            pdf.set_font("Helvetica", "", 8)
            for l in adv_exp_laws:
                pdf.cell(25, 7, l['law_id'], 1, 0, "C")
                pdf.cell(60, 7, l['title'], 1)
                pdf.cell(55, 7, str(l['observed']), 1, 0, "C")
                pdf.cell(50, 7, "Advisory (Experimental)", 1, 1, "C")
            pdf.ln(5)

        # Heuristic signals
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 8, "HEURISTIC ADVISORY SIGNALS", ln=True)
        pdf.ln(2)

        pdf.set_font("Helvetica", "B", 8)
        pdf.set_fill_color(245, 245, 240)
        for i, h in enumerate(hdrs):
            pdf.cell(wds[i], 8, h, 1, 0, "C", True)
        pdf.ln(8)
        pdf.set_font("Helvetica", "", 7)
        for l in heur_laws:
            pdf.cell(18, 7, l['law_id'], 1, 0, "C")
            pdf.cell(52, 7, l['title'], 1)
            pdf.cell(35, 7, str(l['observed']), 1, 0, "C")
            pdf.cell(35, 7, f"{l['operator']} {l['threshold']}", 1, 0, "C")
            pdf.cell(25, 7, l['deviation'], 1, 0, "C")
            pdf.cell(25, 7, pdf.safe(l['units']), 1, 1, "C")

        # ══════ PAGE 4: STRUCTURAL CHARACTERIZATION ══════
        pdf.add_page()
        pdf.section_title("Structural Characterization")

        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, f"Source Type: {source_type.upper()}", ln=True)
        if ss.get("resolution"):
            pdf.cell(w, 8, f"Resolution: {ss['resolution']} A", ln=True)
        pdf.ln(3)

        pdf.cell(w, 8, f"Secondary Structure (Method: {ss.get('method')})", ln=True)
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(63, 10, f"Alpha Helix: {ss.get('helix')}%", 1, 0, "C")
        pdf.cell(63, 10, f"Beta Sheet: {ss.get('sheet')}%", 1, 0, "C")
        pdf.cell(64, 10, f"Loop/Coil: {ss.get('loop')}%", 1, 1, "C")
        pdf.set_font("Helvetica", "I", 7)
        pdf.multi_cell(w, 5, f"Analysis computed over core residues (n={ss.get('core_evaluated')}).")
        pdf.ln(5)

        pdf.set_font("Helvetica", "B", 10)
        dist = conf_meta.get("distribution", {})
        if conf_meta.get("source_type") == "experimental_bfactor":
            pdf.cell(w, 8, "B-factor Distribution", ln=True)
            if dist:
                pdf.set_font("Helvetica", "", 9)
                pdf.cell(63, 10, f"Low (<20): {dist.get('low_bfactor_pct', 'N/A')}%", 1, 0, "C")
                pdf.cell(63, 10, f"Med (20-50): {dist.get('med_bfactor_pct', 'N/A')}%", 1, 0, "C")
                pdf.cell(64, 10, f"High (>50): {dist.get('high_bfactor_pct', 'N/A')}%", 1, 1, "C")
        else:
            pdf.cell(w, 8, "Confidence Distribution (pLDDT)", ln=True)
            if dist:
                pdf.set_font("Helvetica", "", 9)
                pdf.cell(63, 10, f"Core (>=70): {dist.get('high_conf_pct', 'N/A')}%", 1, 0, "C")
                pdf.cell(63, 10, f"Boundary (50-70): {dist.get('med_conf_pct', 'N/A')}%", 1, 0, "C")
                pdf.cell(64, 10, f"Fringe (<50): {dist.get('low_conf_pct', 'N/A')}%", 1, 1, "C")

        # ══════ PAGE 5: PRIORITIZATION MODEL ══════
        pdf.add_page()
        pdf.section_title("Prioritization Model Expansion")
        pdf.set_font("Courier", "B", 10)
        pdf.cell(w, 10, f"FORMULA: {BAYESIAN_FORMULA}", ln=True, align="C")
        pdf.ln(5)
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "Numeric Expansion:", ln=True)
        pdf.set_font("Courier", "", 10)
        pdf.cell(w, 8, f"S6 (Integrity Factor) = {math_data.get('s6', 0.0)}", ln=True)
        pdf.cell(w, 8, f"W_arch ({math_data.get('architecture')}) = {math_data.get('w_arch', 1.0)}", ln=True)
        pdf.cell(w, 8, f"M_S8 (Negative Penalty) = {math_data.get('m_s8', 0.0)}", ln=True)
        pdf.ln(5)
        if v.get("suppression_reason"):
            pdf.set_font("Helvetica", "B", 10)
            pdf.set_text_color(170, 35, 35)
            pdf.cell(w, 10, f"OVERRIDE: {v['suppression_reason']}", border=1, ln=True, align="C")
            pdf.set_text_color(0, 0, 0)
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(w, 15, f"COMPOSITE SCORE: {payload.get('tier3', {}).get('probability', 0)}%", border=1, align="C")

        # ══════ PAGE 6: PROVENANCE ══════
        pdf.add_page()
        pdf.section_title("Provenance & Authentication")
        pdf.set_font("Courier", "B", 9)
        pdf.cell(w, 8, f"AUDIT_ID:       {gov.get('audit_id')}", ln=True)
        pdf.cell(w, 8, f"SOURCE:         {prov.get('source')}", ln=True)
        pdf.cell(w, 8, f"SOURCE_TYPE:    {source_type}", ln=True)
        pdf.cell(w, 8, f"TIMESTAMP_UTC:  {gov.get('timestamp_utc')}", ln=True)
        pdf.cell(w, 8, f"TOTAL_RESIDUES: {ss.get('total_residues')}", ln=True)
        pdf.cell(w, 8, f"TOTAL_ATOMS:    {ss.get('total_atoms')}", ln=True)
        if ss.get("resolution"):
            pdf.cell(w, 8, f"RESOLUTION:     {ss['resolution']} A", ln=True)
        pdf.ln(5)
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "Input Authentication (SHA-256):", ln=True)
        pdf.set_font("Courier", "", 8)
        pdf.multi_cell(w, 4, prov.get('hash', 'N/A'))
        pdf.ln(5)
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(w, 8, "Canon Fingerprint:", ln=True)
        pdf.set_font("Courier", "", 8)
        pdf.cell(w, 8, LAW_CANON_HASH, ln=True)

        pdf.ln(10)
        pdf.set_font("Helvetica", "B", 9)
        pdf.set_fill_color(240, 240, 240)
        if has_advisory:
            pdf.multi_cell(w, 5,
                "PIL-CAL-02: LAW-100 (Bond Integrity) is classified as Advisory for experimental structures. "
                "Engh-Huber sigma values describe restraint targets, not coordinate uncertainty at "
                "macromolecular resolution. RMSZ is reported for transparency but excluded from "
                "deterministic adjudication. Experimental structures have passed wwPDB validation.",
                border=1, fill=True
            )

        return force_bytes(pdf.output(dest='S'))
    except Exception as e:
        err = FPDF()
        err.add_page()
        err.set_font("Helvetica", "B", 12)
        err.cell(0, 10, f"CRITICAL REPORTING ERROR: {str(e)}", ln=True)
        return force_bytes(err.output(dest='S'))
