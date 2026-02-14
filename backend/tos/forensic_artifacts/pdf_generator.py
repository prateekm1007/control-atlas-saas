from fpdf import FPDF
from datetime import datetime
from ..utils.type_guards import sanitize_notary_text, force_bytes
from ..governance.station_sop import (
    LAW_CANON, SCORE_DEFINITIONS, STATION_METADATA, BAYESIAN_FORMULA,
    LAW_CATEGORIES, CATEGORY_DISPLAY, get_laws_by_category, LAW_CANON_HASH
)

class ToscaniniDossier(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 7)
        self.set_text_color(120, 120, 120)
        self.cell(95, 8, f"TOSCANINI v{STATION_METADATA['version']}", align="L")
        self.cell(0, 8, datetime.now().strftime('%Y-%m-%d %H:%M UTC'), align="R", ln=True)
        self.set_draw_color(180, 180, 180)
        self.line(10, 16, 200, 16)
        self.ln(4)
    def footer(self):
        self.set_y(-12)
        self.set_font("Helvetica", "I", 6)
        self.set_text_color(150, 150, 150)
        self.cell(0, 8, f"{STATION_METADATA['institution']} // Page {self.page_no()}", align="C")
    def section_title(self, text, size=14):
        self.set_font("Helvetica", "B", size)
        self.set_text_color(0, 0, 0)
        self.cell(0, 10, text, ln=True)
        self.ln(2)
    def safe(self, text):
        return sanitize_notary_text(str(text))

def generate_v21_dossier(payload):
    try:
        pdf = ToscaniniDossier()
        pdf.set_auto_page_break(True, margin=15)
        w = 190
        v = payload.get('verdict', {})
        ss = payload.get("characterization", {})
        t3 = payload.get("tier3", {})
        prov = payload.get("provenance", {})
        conf_meta = payload.get("confidence_meta", {})
        nkg_data = payload.get("nkg_assessment", {})
        bcomp = payload.get("bayesian_components", {})
        ai_model = payload.get("ai_model_used", "Not specified")
        binary = v.get('binary', 'ERROR')
        phys_score = v.get('physical_score', 0)
        det_score = v.get('deterministic_score', phys_score)
        adv_score = v.get('advisory_score', 100)
        conf_score = v.get('confidence_score', 0)
        conf_available = v.get('confidence_available', True)
        epi = t3.get('probability', 0)
        laws_passed = v.get('laws_passed', 0)
        laws_total = v.get('laws_total', len(LAW_CANON))
        det_passed = v.get('det_passed', laws_passed)
        det_total = v.get('det_total', laws_total)
        heur_passed = v.get('heur_passed', 0)
        heur_total = v.get('heur_total', 0)

        # ══════ PAGE 1: AUDIT REPORT ══════
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 18)
        pdf.cell(w, 14, "STRUCTURAL AUDIT REPORT", ln=True, align="C")
        pdf.ln(2)
        pdf.set_font("Helvetica", "", 8)
        pdf.set_text_color(60, 60, 60)
        pdf.cell(w, 5, f"Program: {payload.get('governance',{}).get('program_id','N/A')}    Artifact: {prov.get('audit_id','N/A')}    Source: {pdf.safe(prov.get('source','N/A'))}", ln=True)
        pdf.set_text_color(0, 0, 0)
        pdf.ln(4)

        # Verdict banner
        if binary == "PASS":
            pdf.set_fill_color(15, 110, 55)
        elif binary == "VETO":
            pdf.set_fill_color(170, 35, 35)
        else:
            pdf.set_fill_color(100, 100, 100)
        pdf.set_text_color(255, 255, 255)
        pdf.set_font("Helvetica", "B", 16)
        pdf.cell(w, 13, f"VERDICT: {binary}", ln=True, align="C", fill=True)
        pdf.set_text_color(0, 0, 0)
        pdf.ln(6)

        # 4-metric grid
        c1w, c2w, c3w, c4w = 50, 42, 50, 48
        pdf.set_font("Helvetica", "B", 7)
        pdf.set_text_color(100, 100, 100)
        pdf.cell(c1w, 5, "DETERMINISTIC", 0, 0, "C")
        pdf.cell(c2w, 5, "ADVISORY", 0, 0, "C")
        conf_label = "pLDDT" if "plddt" in str(conf_meta.get('source_type','')).lower() else "CONFIDENCE"
        pdf.cell(c3w, 5, conf_label, 0, 0, "C")
        pdf.cell(c4w, 5, "EPI PRIORITY", 0, 1, "C")
        pdf.set_text_color(0, 0, 0)
        pdf.set_font("Helvetica", "B", 20)
        pdf.cell(c1w, 12, f"{det_score}%", 0, 0, "C")
        pdf.cell(c2w, 12, f"{adv_score}%", 0, 0, "C")
        if conf_available:
            pdf.cell(c3w, 12, f"{conf_score}", 0, 0, "C")
        else:
            pdf.set_text_color(140, 110, 30)
            pdf.cell(c3w, 12, "N/A", 0, 0, "C")
            pdf.set_text_color(0, 0, 0)
        pdf.cell(c4w, 12, f"{epi}%", 0, 1, "C")
        pdf.set_font("Helvetica", "", 6)
        pdf.set_text_color(100, 100, 100)
        pdf.cell(c1w, 4, f"{det_passed}/{det_total} laws (VETO gate)", 0, 0, "C")
        pdf.cell(c2w, 4, f"{heur_passed}/{heur_total} (informational)", 0, 0, "C")
        pdf.cell(c3w, 4, "Excluded from EPI" if not conf_available else pdf.safe(conf_meta.get('provenance_method','auto')), 0, 0, "C")
        pdf.cell(c4w, 4, "Heuristic composite", 0, 1, "C")
        pdf.set_text_color(0, 0, 0)
        pdf.ln(4)

        # Score definitions
        for key in ["DETERMINISTIC_SCORE", "ADVISORY_SCORE", "ML_CONFIDENCE", "STRATEGIC_SCORE"]:
            d = SCORE_DEFINITIONS.get(key, {})
            if not d: continue
            pdf.set_font("Helvetica", "B", 8)
            pdf.cell(w, 4, d.get('title', key), ln=True)
            pdf.set_font("Helvetica", "", 7)
            pdf.multi_cell(w, 3.5, pdf.safe(f"{d.get('explanation','')} {d.get('impact','')}"))
            pdf.ln(1)

        # Governance + scope
        pdf.ln(2)
        pdf.set_draw_color(180, 180, 180)
        pdf.line(10, pdf.get_y(), 200, pdf.get_y())
        pdf.ln(2)
        pdf.set_font("Helvetica", "B", 8)
        pdf.cell(w, 5, "VETO AUTHORITY POLICY", ln=True)
        pdf.set_font("Helvetica", "", 7)
        pdf.multi_cell(w, 3.5, f"Only deterministic laws ({det_total} of {laws_total}) can trigger a VETO. Heuristic laws ({heur_total} of {laws_total}) are advisory only.")
        if not conf_available:
            pdf.set_font("Helvetica", "I", 7)
            pdf.set_text_color(140, 110, 30)
            pdf.multi_cell(w, 3.5, "NOTICE: ML confidence metadata unavailable. ML component excluded from EPI composite (not penalized to zero).")
            pdf.set_text_color(0, 0, 0)
        pdf.ln(2)
        pdf.set_font("Helvetica", "B", 8)
        pdf.cell(w, 5, "SCOPE & LIMITATIONS", ln=True)
        pdf.set_font("Helvetica", "", 7)
        pdf.multi_cell(w, 3.5, "This audit evaluates geometric and chemical invariants from atomic coordinates. It does not assess biological function, binding affinity, or therapeutic efficacy.")
        pdf.ln(1)
        pdf.set_font("Helvetica", "B", 8)
        pdf.cell(w, 5, "DETERMINISM DECLARATION", ln=True)
        pdf.set_font("Helvetica", "", 7)
        pdf.multi_cell(w, 3.5, f"All DETERMINISTIC law measurements are reproducible given identical input bytes and Toscanini v{STATION_METADATA['version']}. Heuristic measurements use approximate algorithms. AI narrative sections are non-deterministic.")

        # ══════ PAGE 2: TIER-1 AUDIT ══════
        pdf.add_page()
        pdf.section_title(f"TIER-1 PHYSICAL INVARIANT AUDIT ({laws_total} Laws)")
        pdf.set_font("Helvetica", "", 7)
        pdf.cell(w, 4, f"Deterministic: {det_passed}/{det_total} passed  |  Advisory: {heur_passed}/{heur_total} passed  |  Deterministic threshold: 100% required", ln=True)
        pdf.set_font("Helvetica", "I", 6)
        pdf.set_text_color(100, 100, 100)
        pdf.cell(w, 3, f"Combined physical score (legacy): {phys_score}% ({laws_passed}/{laws_total} total)", ln=True)
        pdf.set_text_color(0, 0, 0)
        pdf.ln(3)

        all_laws_list = payload.get('tier1', {}).get('laws', [])
        laws_by_id = {l.get('law_id'): l for l in all_laws_list}
        cat_groups = get_laws_by_category()

        for cat_key in ["geometric", "chemical"]:
            section_name = CATEGORY_DISPLAY.get(cat_key, cat_key.upper())
            cat_lids = cat_groups.get(cat_key, [])
            section_laws = [laws_by_id[lid] for lid in cat_lids if lid in laws_by_id]
            pdf.set_font("Helvetica", "B", 10)
            pdf.set_text_color(40, 40, 40)
            pdf.cell(w, 7, section_name, ln=True)
            pdf.set_text_color(0, 0, 0)

            for l in section_laws:
                if pdf.get_y() > 260:
                    pdf.add_page()
                lid = l.get('law_id', 'N/A')
                status = l.get('status', 'N/A')
                measurement = l.get('measurement', '')
                method = l.get('method', 'unknown')
                title = l.get('title', lid)
                if status == "PASS":
                    pdf.set_text_color(15, 110, 55)
                elif status == "FAIL":
                    pdf.set_text_color(170, 35, 35)
                else:
                    pdf.set_text_color(100, 100, 100)
                pdf.set_font("Helvetica", "B", 8)
                pdf.cell(16, 5, f"[{status}]", 0, 0)
                pdf.set_text_color(0, 0, 0)
                pdf.cell(0, 5, f"{lid}: {pdf.safe(title)}", ln=True)
                pdf.set_font("Helvetica", "", 7)
                pdf.set_x(pdf.l_margin + 16)
                pdf.multi_cell(w - 16, 3.5, pdf.safe(str(measurement)))
                if method == "heuristic":
                    pdf.set_font("Helvetica", "I", 6)
                    pdf.set_text_color(140, 110, 30)
                    pdf.set_x(pdf.l_margin + 16)
                    pdf.cell(0, 3, "HEURISTIC METHOD - result is approximate", ln=True)
                    pdf.set_text_color(0, 0, 0)
                pdf.ln(1)
            pdf.ln(2)

        pdf.set_font("Helvetica", "I", 6)
        pdf.set_text_color(100, 100, 100)
        pdf.cell(w, 3, "Deterministic = exact geometric computation. Heuristic = approximate algorithm.", ln=True)
        pdf.set_text_color(0, 0, 0)

        # ══════ PAGE 3: CHARACTERIZATION ══════
        pdf.add_page()
        pdf.section_title("STRUCTURAL CHARACTERIZATION")
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 6, "Secondary Structure Distribution", ln=True)
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(63, 7, f"Alpha Helix: {ss.get('helix',0)}%", border=1, align="C")
        pdf.cell(63, 7, f"Beta Sheet: {ss.get('sheet',0)}%", border=1, align="C")
        pdf.cell(63, 7, f"Coil/Loop: {ss.get('loop',0)}%", border=1, ln=True, align="C")
        pdf.set_font("Helvetica", "I", 6)
        pdf.set_text_color(100, 100, 100)
        pdf.cell(w, 4, f"Method: {ss.get('method','heuristic')} // Characterization only - not used in VETO or EPI", ln=True, align="C")
        pdf.set_text_color(0, 0, 0)
        pdf.ln(4)

        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 6, "Structural Metrics", ln=True)
        for label, val in [("Total Atoms", ss.get('total_atoms', 'N/A')), ("Total Residues", ss.get('total_residues', 'N/A')), ("Source", pdf.safe(prov.get('source', 'N/A')))]:
            pdf.set_font("Helvetica", "", 8)
            pdf.set_text_color(80, 80, 80)
            pdf.cell(65, 6, str(label), border="B")
            pdf.set_text_color(0, 0, 0)
            pdf.cell(0, 6, str(val), border="B", ln=True)
        pdf.ln(4)

        # Confidence decomposition — aligned with Page 1
        conf_src_lower = str(conf_meta.get('source_type','')).lower()
        conf_unit = "pLDDT (0-100)" if 'plddt' in conf_src_lower else "B-factor (A^2)" if 'b-factor' in conf_src_lower else "(unitless)"
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 6, "Confidence Decomposition", ln=True)
        pdf.set_font("Helvetica", "", 7)
        pdf.cell(w, 4, f"Source: {pdf.safe(conf_meta.get('source_type','unknown'))}  |  Detection: {conf_meta.get('provenance_method','auto')}", ln=True)
        pdf.ln(2)

        conf_data_avail = conf_meta.get('data_available', False)
        if conf_data_avail and conf_meta.get('mean', 0) > 0:
            for label, val in [(f"Mean {conf_unit}", conf_meta.get('mean')), (f"Min {conf_unit}", conf_meta.get('min')), (f"Max {conf_unit}", conf_meta.get('max')), (f"Std Dev {conf_unit}", conf_meta.get('std')), ("Low-confidence (<50)", len(conf_meta.get('low_confidence_regions',[])))]:
                pdf.set_font("Helvetica", "", 8)
                pdf.set_text_color(80, 80, 80)
                pdf.cell(65, 5, str(label), border="B")
                pdf.set_text_color(0, 0, 0)
                pdf.cell(0, 5, str(val), border="B", ln=True)
        else:
            pdf.set_font("Helvetica", "", 8)
            pdf.set_text_color(80, 80, 80)
            pdf.cell(65, 5, f"Mean {conf_unit}", border="B")
            pdf.set_text_color(0, 0, 0)
            pdf.cell(0, 5, str(conf_score), border="B", ln=True)
            pdf.set_text_color(80, 80, 80)
            pdf.cell(65, 5, "Per-residue breakdown", border="B")
            pdf.set_text_color(0, 0, 0)
            pdf.cell(0, 5, "Not available for this structure", border="B", ln=True)
            pdf.ln(2)
            pdf.set_font("Helvetica", "I", 7)
            pdf.set_text_color(140, 110, 30)
            pdf.multi_cell(w, 3.5, f"Structure-level mean ({conf_score}) derived from metadata. Used consistently across all pages.")
            pdf.set_text_color(0, 0, 0)

        # ══════ PAGE 4: PRIORITIZATION ══════
        pdf.add_page()
        pdf.section_title("HEURISTIC PRIORITIZATION MODEL (EPI)")
        pdf.set_font("Courier", "", 8)
        pdf.set_fill_color(240, 240, 240)
        pdf.cell(w, 8, f"Formula: {BAYESIAN_FORMULA}", ln=True, align="C", fill=True)
        pdf.ln(4)
        breakdown = bcomp.get("breakdown", [])
        pdf.set_font("Helvetica", "B", 8)
        pdf.cell(8, 6, "#", border="B"); pdf.cell(62, 6, "Component", border="B"); pdf.cell(30, 6, "Value", border="B"); pdf.cell(0, 6, "Derivation", border="B", ln=True)
        for i, comp in enumerate(breakdown, 1):
            pdf.set_font("Helvetica", "", 8)
            pdf.cell(8, 6, str(i), border="B"); pdf.cell(62, 6, str(comp.get('name','')), border="B"); pdf.cell(30, 6, str(comp.get('value',0)), border="B"); pdf.cell(0, 6, pdf.safe(str(comp.get('source',''))), border="B", ln=True)
        pdf.ln(4)
        s6c = next((c['value'] for c in breakdown if 'composite' in str(c.get('name','')).lower()), 0)
        wa = next((c['value'] for c in breakdown if 'w_arch' in str(c.get('name','')).lower()), 1.0)
        ms = next((c['value'] for c in breakdown if 's8' in str(c.get('name','')).lower() or 'nkg' in str(c.get('name','')).lower()), 0.0)
        pdf.set_font("Courier", "B", 9)
        pdf.set_fill_color(240, 240, 240)
        pdf.cell(w, 8, f"P = {s6c} x {wa} x (1 - {ms}) = {round(s6c * wa * (1 - ms), 4)}", ln=True, align="C", fill=True)
        pdf.ln(2)
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(w, 10, f"FINAL EPI INDEX: {epi}%", ln=True, align="C")
        pdf.ln(2)
        pdf.set_font("Helvetica", "I", 7)
        pdf.set_text_color(100, 100, 100)
        pdf.multi_cell(w, 3.5, "S6_phys uses the DETERMINISTIC integrity score. When confidence unavailable, S6_ml is excluded. This is a heuristic composite, not a Bayesian posterior.")
        pdf.set_text_color(0, 0, 0)

        # ══════ PAGE 5: NKG ══════
        pdf.add_page()
        pdf.section_title("NEGATIVE KNOWLEDGE GRAPH ASSESSMENT")
        nkg_records = nkg_data.get('records_searched', 0)
        for label, val in [("Similarity Threshold", nkg_data.get('threshold','N/A')), ("Records Searched", nkg_records), ("Matches Above Threshold", nkg_data.get('matches_above',0)), ("Closest Failure", nkg_data.get('closest_similarity','N/A')), ("S8 Penalty", nkg_data.get('penalty',0.0))]:
            pdf.set_font("Helvetica", "", 8)
            pdf.set_text_color(80, 80, 80)
            pdf.cell(65, 5, str(label))
            pdf.set_text_color(0, 0, 0)
            pdf.cell(0, 5, str(val), ln=True)
        pdf.ln(3)
        pdf.set_font("Helvetica", "", 8)
        if nkg_records == 0:
            pdf.multi_cell(w, 4, "NKG deployment state: initialized. Database contains 0 archived failure records. S8 penalty inactive until sufficient historical data present.")
        else:
            pdf.multi_cell(w, 4, pdf.safe(nkg_data.get('narrative', 'Assessment complete.')))

        # ══════ PAGE 6: PROVENANCE ══════
        pdf.add_page()
        pdf.section_title("DATA PROVENANCE & INTEGRITY")
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 6, "INPUT PROVENANCE", ln=True)
        byte_count = prov.get('byte_count', 0)
        byte_exact = prov.get('byte_count_exact', False)
        for label, val in [("Source", pdf.safe(prov.get('source','N/A'))), ("Format", "mmCIF" if str(prov.get('source','')).lower().endswith(('.cif','.mmcif')) else "PDB"), ("Input Size (bytes)", f"{byte_count:,} ({'exact' if byte_exact else 'estimated'})"), ("Audit ID", prov.get('audit_id','N/A')), ("Station Version", prov.get('station_version', STATION_METADATA['version'])), ("Timestamp", datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC'))]:
            pdf.set_font("Helvetica", "", 8)
            pdf.set_text_color(80, 80, 80)
            pdf.cell(65, 5, str(label))
            pdf.set_text_color(0, 0, 0)
            pdf.cell(0, 5, str(val), ln=True)
        pdf.ln(6)
        pdf.set_font("Courier", "B", 10)
        pdf.cell(w, 8, "SHA-256 INTEGRITY CHECKSUM", ln=True, align="C")
        pdf.set_font("Courier", "", 7)
        pdf.multi_cell(w, 4, prov.get('hash', 'UNAVAILABLE'), align="C")
        pdf.ln(3)
        pdf.set_font("Helvetica", "I", 7)
        pdf.set_text_color(100, 100, 100)
        pdf.multi_cell(w, 3.5, "Checksum computed from raw input bytes. Any file modification produces different checksum.", align="C")
        pdf.set_text_color(0, 0, 0)
        pdf.ln(4)

        # AI disclosure
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 6, "AI-GENERATED CONTENT DISCLOSURE", ln=True, align="C")
        pdf.set_font("Helvetica", "", 7)
        ai_display = pdf.safe(ai_model)
        if "exhausted" in ai_display.lower() or "Static" in ai_display:
            ai_display = "Internal analysis module (deterministic fallback)"
        pdf.multi_cell(w, 3.5, f"Narrative by: {ai_display}. AI text is non-deterministic. All measurements are computed deterministically.", align="C")
        pdf.ln(4)

        # Version immutability
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 6, "VERSION IMMUTABILITY", ln=True, align="C")
        pdf.set_font("Helvetica", "", 7)
        pdf.multi_cell(w, 3.5, f"Results tied to v{STATION_METADATA['version']} ({len(LAW_CANON)} laws). Different version may produce different results.", align="C")
        pdf.ln(4)

        # Canon hash
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(w, 6, "LAW CANON FINGERPRINT", ln=True, align="C")
        pdf.set_font("Courier", "", 7)
        pdf.cell(w, 5, f"Canon Hash: {LAW_CANON_HASH}", ln=True, align="C")
        pdf.set_font("Helvetica", "I", 6)
        pdf.set_text_color(100, 100, 100)
        pdf.cell(w, 4, "SHA-256 over sorted law definitions. Changes if any law is modified.", ln=True, align="C")
        pdf.set_text_color(0, 0, 0)

        return force_bytes(pdf.output(dest='S'))
    except Exception as e:
        err = FPDF()
        err.add_page()
        err.set_font("Helvetica", "B", 28)
        err.set_text_color(200, 30, 30)
        err.cell(0, 30, "GENERATION FAILED", ln=True, align="C")
        err.set_font("Helvetica", "", 9)
        err.set_text_color(0, 0, 0)
        err.multi_cell(0, 5, f"Error: {str(e)[:500]}")
        return force_bytes(err.output(dest='S'))
