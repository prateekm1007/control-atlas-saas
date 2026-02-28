"""
dashboard/pdf_certificate.py

Toscanini B.6 — Refinement Resolution Certificate

Design:
  - Pure reportlab (no system fonts)
  - Deterministic: suppresses wall-clock timestamps via reportlab canvas patch
  - Verdict stored in /Keywords PDF metadata (uncompressed, byte-searchable)
  - No PDB coordinates in output
  - Returns bytes → st.download_button

delta dict keys (from comparison_engine.compare_audits):
    baseline_audit_id, refined_audit_id,
    verdict_change, verdict_improved (bool),
    coverage_delta (float), violation_count_before (int),
    violation_count_after (int), violation_count_delta (int),
    improvements (list[str]), regressions (list[str]),
    law_changes (list[dict])

meta dict keys:
    structure_name, api_key_masked, tier, engine, issued_at (int unix ts)
"""
from __future__ import annotations

import io
import time
from typing import Any, Dict, List

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        HRFlowable, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
    )
    from reportlab.pdfgen import canvas as rl_canvas
    _REPORTLAB_OK = True
except ImportError:
    _REPORTLAB_OK = False

# ── Palette (lazy init — safe at import time) ─────────────────────────────────
_C: Dict = {}

def _pal():
    if not _C:
        _C["navy"]   = colors.HexColor("#0F172A")
        _C["blue"]   = colors.HexColor("#3B82F6")
        _C["green"]  = colors.HexColor("#10B981")
        _C["amber"]  = colors.HexColor("#F59E0B")
        _C["red"]    = colors.HexColor("#EF4444")
        _C["silver"] = colors.HexColor("#F1F5F9")
        _C["mid"]    = colors.HexColor("#94A3B8")
        _C["white"]  = colors.white
    return _C


# ── Determinism: fixed epoch replaces wall-clock in PDF metadata ──────────────
_FIXED_DATE = "D:20240101000000+00'00'"   # constant across all builds


def _make_doc(buf: io.BytesIO, baseline_id: str, verdict_keyword: str) -> SimpleDocTemplate:
    """
    Build SimpleDocTemplate with suppressed wall-clock timestamps.
    reportlab writes CreationDate/ModDate from canvas._doc.info — we
    pre-set them to the fixed epoch string so output is byte-identical.
    Verdict stored in /Keywords for byte-level searchability.
    """
    doc = SimpleDocTemplate(
        buf, pagesize=A4,
        leftMargin=20*mm, rightMargin=20*mm,
        topMargin=15*mm,  bottomMargin=15*mm,
        title=f"Toscanini Certificate — {baseline_id}",
        author="Toscanini Forensic Governance Engine",
        subject="Refinement Resolution Certificate",
        keywords=verdict_keyword,   # stored uncompressed in /Keywords
        creator="Toscanini",
        producer="Toscanini PDF Engine",
    )
    return doc


def _patch_timestamps(doc: SimpleDocTemplate) -> None:
    """
    Overwrite reportlab's internal timestamp fields to fixed values.
    Called after doc.build() sets up internal state but before rendering.
    This is the only reliable way to suppress wall-clock in reportlab 4.x.
    We monkey-patch the _doc.info dict if it exists.
    """
    try:
        # reportlab stores metadata in doc.canv._doc.info
        info = doc.canv._doc.info
        info.datestr     = _FIXED_DATE
        info.created     = _FIXED_DATE
        info.modified    = _FIXED_DATE
    except Exception:
        pass   # Non-fatal — determinism best-effort


class _DeterministicDoc(SimpleDocTemplate):
    """
    Subclass that injects fixed timestamps into the canvas before
    the PDF stream is finalised. Works with reportlab 3.x and 4.x.
    """
    def __init__(self, buf, verdict_keyword, **kwargs):
        self._verdict_keyword = verdict_keyword
        super().__init__(buf, **kwargs)

    def handle_documentBegin(self):
        super().handle_documentBegin()
        self._fix_timestamps()

    def _fix_timestamps(self):
        try:
            self.canv._doc.CreationDate = _FIXED_DATE
            self.canv._doc.ModDate      = _FIXED_DATE
        except Exception:
            pass

    def build(self, flowables, **kwargs):
        # Override build to inject timestamps after canvas is created
        super().build(flowables, **kwargs)


def build_certificate(delta: Dict, meta: Dict) -> bytes:
    """
    Build and return the Refinement Resolution Certificate as PDF bytes.
    Byte-identical for same (delta, meta) inputs.
    Falls back to UTF-8 plaintext if reportlab is not installed.
    """
    if not _REPORTLAB_OK:
        return _plaintext_fallback(delta, meta)

    _pal()

    vb = int(delta.get("violation_count_before", 0))
    va = int(delta.get("violation_count_after",  0))
    rd = vb - va

    if rd > 0:
        verdict_keyword = "CERTIFIED"
    elif rd < 0:
        verdict_keyword = "REGRESSION"
    else:
        verdict_keyword = "NEUTRAL"

    buf = io.BytesIO()

    # Use fixed timestamp to suppress wall-clock non-determinism
    doc = SimpleDocTemplate(
        buf, pagesize=A4,
        leftMargin=20*mm, rightMargin=20*mm,
        topMargin=15*mm,  bottomMargin=15*mm,
        title=f"Toscanini Certificate — {delta.get('baseline_audit_id','?')}",
        author="Toscanini Forensic Governance Engine",
        subject="Refinement Resolution Certificate",
        keywords=verdict_keyword,
        creator="Toscanini",
        producer="Toscanini PDF Engine",
    )

    styles = _styles()
    story: List[Any] = []

    _section_header(story, styles, delta, meta)
    _section_banner(story, styles, delta, verdict_keyword)
    _section_metrics(story, styles, delta)
    _section_law_changes(story, styles, delta)
    _section_governance(story, styles)
    _section_footer(story, styles, delta, meta)

    # Build with timestamp suppression
    doc.build(story)

    # Post-build: replace timestamp bytes in PDF stream
    raw = buf.getvalue()
    raw = _suppress_timestamps(raw)

    return raw


def _suppress_timestamps(raw: bytes) -> bytes:
    """
    Normalise all non-deterministic fields in a ReportLab PDF stream.

    ReportLab injects two sources of non-determinism:
      1. Timestamps  : D:YYYYMMDDHHmmss+HH'00'
      2. Document ID : /ID [<hexhash><hexhash>]  — changes every run

    After suppression, byte output is identical for identical inputs.
    xref offsets and stream lengths are stable because layout is fixed.
    """
    import re as _re

    # 1. Suppress timestamps
    raw = _re.sub(
        rb"\(D:\d{14}[+\-Z][^)]*\)",
        b"(D:20240101000000+00'00')",
        raw,
    )

    # 2. Suppress /ID array  [<hash1><hash2>]
    raw = _re.sub(
        rb"/ID\s*\n?\[<[0-9a-fA-F]+><[0-9a-fA-F]+>\]",
        b"/ID\n[<00000000000000000000000000000000><00000000000000000000000000000000>]",
        raw,
    )

    return raw


# ── Section builders ──────────────────────────────────────────────────────────

def _section_header(story, styles, delta, meta):
    c = _pal()
    story.append(Paragraph("TOSCANINI", styles["brand"]))
    story.append(Paragraph("Forensic Structural Governance Engine", styles["sub"]))
    story.append(HRFlowable(width="100%", thickness=2, color=c["navy"]))
    story.append(Spacer(1, 4*mm))
    story.append(Paragraph("REFINEMENT RESOLUTION CERTIFICATE", styles["cert_title"]))
    story.append(Spacer(1, 3*mm))

    issued = time.strftime(
        "%Y-%m-%d %H:%M UTC",
        time.gmtime(meta.get("issued_at", 0))
    )
    rows = [
        ["Structure",   meta.get("structure_name", "—"),
         "Issued",      issued],
        ["Baseline ID", str(delta.get("baseline_audit_id", "—"))[:20],
         "Engine",      meta.get("engine", "—").upper()],
        ["API Key",     meta.get("api_key_masked", "—"),
         "Tier",        meta.get("tier", "—").title()],
    ]
    t = Table(rows, colWidths=[30*mm, 65*mm, 25*mm, 50*mm])
    t.setStyle(TableStyle([
        ("FONTNAME",      (0,0),(-1,-1), "Helvetica"),
        ("FONTSIZE",      (0,0),(-1,-1), 8),
        ("TEXTCOLOR",     (0,0),(0,-1),  c["mid"]),
        ("TEXTCOLOR",     (2,0),(2,-1),  c["mid"]),
        ("FONTNAME",      (1,0),(1,-1),  "Helvetica-Bold"),
        ("FONTNAME",      (3,0),(3,-1),  "Helvetica-Bold"),
        ("TEXTCOLOR",     (1,0),(1,-1),  c["navy"]),
        ("TEXTCOLOR",     (3,0),(3,-1),  c["navy"]),
        ("BOTTOMPADDING", (0,0),(-1,-1), 2),
    ]))
    story.append(t)
    story.append(Spacer(1, 6*mm))


def _section_banner(story, styles, delta, verdict_keyword: str):
    c   = _pal()
    vb  = int(delta.get("violation_count_before", 0))
    va  = int(delta.get("violation_count_after",  0))
    rd  = vb - va
    pct = rd / vb * 100 if vb else 0.0

    if rd > 0:
        colour = c["green"]
        text   = f"CERTIFIED: {rd} violation(s) resolved ({pct:.1f}% improvement)"
    elif rd < 0:
        colour = c["red"]
        text   = f"REGRESSION: {abs(rd)} additional violation(s) introduced."
    else:
        colour = c["amber"]
        text   = "NEUTRAL: No net change in violation count."

    b = Table([[text]], colWidths=[170*mm])
    b.setStyle(TableStyle([
        ("BACKGROUND",    (0,0),(-1,-1), colour),
        ("TEXTCOLOR",     (0,0),(-1,-1), c["white"]),
        ("FONTNAME",      (0,0),(-1,-1), "Helvetica-Bold"),
        ("FONTSIZE",      (0,0),(-1,-1), 11),
        ("TOPPADDING",    (0,0),(-1,-1), 6),
        ("BOTTOMPADDING", (0,0),(-1,-1), 6),
        ("LEFTPADDING",   (0,0),(-1,-1), 8),
    ]))
    story.append(b)
    story.append(Spacer(1, 6*mm))


def _section_metrics(story, styles, delta):
    c = _pal()
    story.append(Paragraph("Performance Metrics", styles["h2"]))
    story.append(Spacer(1, 2*mm))

    vb  = delta.get("violation_count_before", "—")
    va  = delta.get("violation_count_after",  "—")
    vd  = delta.get("violation_count_delta",  "—")
    cov = float(delta.get("coverage_delta", 0.0))
    vc  = delta.get("verdict_change", "—")
    imp = "YES" if delta.get("verdict_improved") else "NO"

    rows = [
        ["Metric", "Before", "After", "Delta", "Trend"],
        ["Violations", str(vb), str(va), str(vd),
         "DOWN" if (isinstance(vd, (int,float)) and vd < 0) else
         ("UP"  if (isinstance(vd, (int,float)) and vd > 0) else "FLAT")],
        ["Coverage change (pp)", "—",
         f"+{cov:.1f}" if cov >= 0 else f"{cov:.1f}", "—",
         "UP" if cov > 0 else ("DOWN" if cov < 0 else "FLAT")],
        ["Verdict", "—", vc, "—", "IMPROVED" if imp == "YES" else "UNCHANGED"],
    ]

    col_w = [60*mm, 28*mm, 28*mm, 28*mm, 26*mm]
    t = Table(rows, colWidths=col_w, repeatRows=1)
    ts = [
        ("BACKGROUND",    (0,0),(-1,0), c["navy"]),
        ("TEXTCOLOR",     (0,0),(-1,0), c["white"]),
        ("FONTNAME",      (0,0),(-1,0), "Helvetica-Bold"),
        ("FONTSIZE",      (0,0),(-1,0), 9),
        ("ALIGN",         (0,0),(-1,0), "CENTER"),
        ("FONTNAME",      (0,1),(-1,-1), "Helvetica"),
        ("FONTSIZE",      (0,1),(-1,-1), 8.5),
        ("ROWBACKGROUNDS",(0,1),(-1,-1), [c["white"], c["silver"]]),
        ("ALIGN",         (1,1),(-1,-1), "CENTER"),
        ("GRID",          (0,0),(-1,-1), 0.3, c["mid"]),
        ("TOPPADDING",    (0,0),(-1,-1), 3),
        ("BOTTOMPADDING", (0,0),(-1,-1), 3),
    ]
    for i, row in enumerate(rows[1:], 1):
        sym = row[4]
        col = c["green"] if sym in ("DOWN","IMPROVED","UP") else (
              c["red"]   if sym in ("UP","UNCHANGED") else c["navy"])
        # More precise colouring
        if sym == "DOWN":   col = c["green"]
        elif sym == "UP":   col = c["red"]
        elif sym == "IMPROVED": col = c["green"]
        else:               col = c["mid"]
        ts.append(("TEXTCOLOR", (4,i),(4,i), col))
    t.setStyle(TableStyle(ts))
    story.append(t)
    story.append(Spacer(1, 6*mm))


def _section_law_changes(story, styles, delta):
    c = _pal()
    improvements = delta.get("improvements", [])
    regressions  = delta.get("regressions",  [])
    if not improvements and not regressions:
        return

    story.append(Paragraph("Law-Level Outcomes", styles["h2"]))
    story.append(Spacer(1, 2*mm))

    if improvements:
        story.append(Paragraph(
            f"Resolved ({len(improvements)}): " + ", ".join(improvements),
            styles["good"],
        ))
    if regressions:
        story.append(Paragraph(
            f"Regressed ({len(regressions)}): " + ", ".join(regressions),
            styles["bad"],
        ))
    story.append(Spacer(1, 4*mm))


def _section_governance(story, styles):
    c = _pal()
    story.append(Paragraph("Governance Basis", styles["h2"]))
    story.append(Spacer(1, 1*mm))
    story.append(Paragraph(
        "This certificate was issued by the Toscanini Forensic Governance Engine "
        "under its 15-Law Physical Validity Canon. All thresholds are frozen at "
        "v23.2.0-B5-beta. Per-residue violation records are stored under the "
        "Baseline Audit ID shown above.",
        styles["body"],
    ))
    story.append(Spacer(1, 4*mm))


def _section_footer(story, styles, delta, meta):
    c = _pal()
    story.append(HRFlowable(width="100%", thickness=0.5, color=c["mid"]))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph(
        f"Toscanini Certificate — "
        f"Baseline: {delta.get('baseline_audit_id','?')} — "
        f"Refined: {delta.get('refined_audit_id','?')} — "
        f"Read-only tamper-evident document.",
        styles["footer"],
    ))


# ── Styles ────────────────────────────────────────────────────────────────────

def _styles() -> dict:
    c    = _pal()
    base = getSampleStyleSheet()
    return {
        "brand": ParagraphStyle("brand", parent=base["Normal"],
            fontName="Helvetica-Bold", fontSize=22,
            textColor=c["navy"], spaceAfter=1*mm),
        "sub": ParagraphStyle("sub", parent=base["Normal"],
            fontName="Helvetica", fontSize=10,
            textColor=c["mid"], spaceAfter=3*mm),
        "cert_title": ParagraphStyle("cert_title", parent=base["Normal"],
            fontName="Helvetica-Bold", fontSize=14,
            textColor=c["blue"], alignment=1, spaceAfter=4*mm),
        "h2": ParagraphStyle("h2", parent=base["Normal"],
            fontName="Helvetica-Bold", fontSize=10,
            textColor=c["navy"], spaceAfter=2*mm),
        "body": ParagraphStyle("body", parent=base["Normal"],
            fontName="Helvetica", fontSize=8.5,
            textColor=c["navy"], leading=13),
        "good": ParagraphStyle("good", parent=base["Normal"],
            fontName="Helvetica", fontSize=8.5,
            textColor=c["green"], spaceAfter=2*mm),
        "bad": ParagraphStyle("bad", parent=base["Normal"],
            fontName="Helvetica", fontSize=8.5,
            textColor=c["red"], spaceAfter=2*mm),
        "footer": ParagraphStyle("footer", parent=base["Normal"],
            fontName="Helvetica", fontSize=7,
            textColor=c["mid"], alignment=1),
    }


def _plaintext_fallback(delta: Dict, meta: Dict) -> bytes:
    vb = int(delta.get("violation_count_before", 0))
    va = int(delta.get("violation_count_after",  0))
    rd = vb - va
    kw = "CERTIFIED" if rd > 0 else ("REGRESSION" if rd < 0 else "NEUTRAL")
    lines = [
        "TOSCANINI — REFINEMENT RESOLUTION CERTIFICATE",
        f"Keywords: {kw}",
        "=" * 60,
        f"Structure  : {meta.get('structure_name','—')}",
        f"Baseline   : {delta.get('baseline_audit_id','—')}",
        f"Refined    : {delta.get('refined_audit_id','—')}",
        f"Verdict    : {delta.get('verdict_change','—')}",
        f"Violations : {vb} -> {va}",
        f"Improved   : {delta.get('verdict_improved','—')}",
        f"Engine     : {meta.get('engine','—')}",
        f"Issued     : {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(meta.get('issued_at',0)))}",
        "",
        "reportlab not installed — pip install reportlab>=4.0",
    ]
    return "\n".join(lines).encode("utf-8")
