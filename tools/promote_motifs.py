import shutil
from pathlib import Path

# === SOVEREIGN PROMOTION TABLE (Week 1) ===
PROMOTIONS = {
    # VERIFIED — cleared geometry + stable energetics
    "pdl1_v1_baseline": ("verified", "LAW-108"),
    "pdl1_v2_02_ywpg": ("verified", "LAW-108"),
    "pdl1_v16_P4_proline_champion": ("verified", "LAW-114, LAW-115"),
    "pdl1_v20_02_alanine_cap": ("verified", "LAW-114"),
    "pdl1_v8_01_ws_core": ("verified", "LAW-108"),

    # KRAS — LAW-105
    "kras_g12d_v5": ("rejected", "LAW-105"),
    "kras_g12d_v6": ("rejected", "LAW-105"),
    "kras_g12d_v7": ("rejected", "LAW-105"),
    "kras_g12d_v8": ("rejected", "LAW-105"),

    # CTLA-4 — LAW-106
    "ctla4_v1_transfer": ("rejected", "LAW-106"),
    "ctla4_v3_01_native": ("rejected", "LAW-106"),
    "ctla4_v3_02_charge": ("rejected", "LAW-106"),
    "ctla4_v3_03_hydrophobic": ("rejected", "LAW-106"),
    "ctla4_v4_cysteine_loop": ("rejected", "LAW-106"),
    "ctla4_v5_12mer_reach": ("rejected", "LAW-106"),
    "ctla4_v6_01_wedge": ("rejected", "LAW-106"),
    "ctla4_v6_02_anchor_transfer": ("rejected", "LAW-106"),

    # PD-L1 globular — LAW-111
    "pdl1_v11_trp_cage_graft": ("rejected", "LAW-111"),
    "pdl1_v12_L0_direct": ("rejected", "LAW-111"),
    "pdl1_v12_L2_short": ("rejected", "LAW-111"),
    "pdl1_v12_L4_medium": ("rejected", "LAW-111"),
    "pdl1_v12_L6_long": ("rejected", "LAW-111"),

    # Aromatic overload — LAW-110
    "pdl1_v5_01_double_trp": ("rejected", "LAW-110"),
    "pdl1_v7_02_double_trp_refinement": ("rejected", "LAW-110"),

    # Proline violation — LAW-108
    "pdl1_v10_01_threonine_pivot": ("rejected", "LAW-108"),
}

def promote():
    base = Path("ledger/motifs")
    pending = base / "pending"

    promoted = 0
    for motif, (status, law) in PROMOTIONS.items():
        src = pending / motif
        dst = base / status / motif

        if not src.exists():
            print(f"⚠️ Skipping missing motif: {motif}")
            continue

        # Update JUDGMENT.md
        j = src / "JUDGMENT.md"
        j.write_text(
            j.read_text().replace(
                "- **Status:** PENDING",
                f"- **Status:** {status.upper()}\n- **Verdict:** {'PASS' if status=='verified' else 'FAIL'}"
            )
        )

        # Overwrite DOCTRINE.md deterministically
        d = src / "DOCTRINE.md"
        d.write_text(
            f"# DOCTRINE — {motif}\n\n"
            f"- **Law Linkage:** {law}\n"
        )

        dst.parent.mkdir(exist_ok=True)
        shutil.move(str(src), str(dst))
        promoted += 1
        print(f"✅ {motif} → {status.upper()} ({law})")

    print(f"\n🏁 Total deterministic promotions: {promoted}")

if __name__ == "__main__":
    promote()
