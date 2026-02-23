"""
Toscanini PDF Theme System
Centralized styling constants for forensic dossier generation.
Single source of truth for all colors, fonts, and layout.
"""

PDF_COLORS = {
    # Verdict colors
    "PASS": (46, 125, 50),           # #2E7D32 - Material Green 700
    "VETO": (211, 47, 47),           # #D32F2F - Material Red 700
    "INDETERMINATE": (117, 117, 117), # #757575 - Material Grey 600
    "ALERT": (245, 124, 0),          # #F57C00 - Material Orange 700

    # Text colors
    "TEXT_PRIMARY": (0, 0, 0),       # Black
    "TEXT_SECONDARY": (80, 80, 80),  # Dark grey
    "TEXT_MUTED": (128, 128, 128),   # Mid grey
    "TEXT_SUBTLE": (100, 100, 100),  # Subtle grey
    "TEXT_WHITE": (255, 255, 255),   # White

    # Banner colors
    "BANNER_VETO": (180, 30, 30),    # Deep red
    "BANNER_ALERT": (200, 150, 30),  # Amber
    "BANNER_PASS": (46, 125, 50),    # Green

    # Background colors
    "BG_LIGHT": (240, 240, 240),     # Light grey
    "BG_LIGHTER": (250, 250, 250),   # Near white
    "LIGHT_GREY": (230, 230, 230),   # Light grey (legacy key)
    "BG_VETO_TABLE": (245, 230, 230), # Light red tint
    "BG_PASS_TABLE": (230, 245, 230), # Light green tint
    "BG_ADV_TABLE": (255, 245, 220),  # Light amber tint
    "BG_SOURCE": (240, 240, 240),    # Source header

    # Section divider colors
    "SECTION_DET": (20, 60, 120),    # Deterministic blue
    "SECTION_ADV": (140, 100, 20),   # Advisory amber
}

PDF_LAYOUT = {
    "MARGIN_LEFT": 15,
    "MARGIN_RIGHT": 15,
    "LINE_HEIGHT": 6,
    "LINE_HEIGHT_WIDE": 8,
    "PAGE_WIDTH": 190,
}
