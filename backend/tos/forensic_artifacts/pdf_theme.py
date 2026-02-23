"""
Toscanini PDF Theme System
Centralized styling constants for forensic dossier generation.
"""

PDF_COLORS = {
    "PASS": (46, 125, 50),        # #2E7D32 - Material Green 700
    "VETO": (211, 47, 47),        # #D32F2F - Material Red 700
    "INSUFFICIENT": (117, 117, 117), # #757575 - Material Grey 600
    "ALERT": (245, 124, 0),       # #F57C00 - Material Orange 700
    "HEADER": (0, 0, 0),
    "LIGHT_GREY": (230, 230, 230),
}

PDF_LAYOUT = {
    "MARGIN_LEFT": 15,
    "MARGIN_RIGHT": 15,
    "LINE_HEIGHT": 6,
    "PAGE_WIDTH": 190,  # A4 width minus margins
}
