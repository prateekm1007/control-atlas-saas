from typing import Any

def force_bytes(data: Any) -> bytes:
    if isinstance(data, (bytes, bytearray)): return bytes(data)
    return str(data).encode('utf-8', errors='ignore')

def force_str(data: Any) -> str:
    if isinstance(data, str): return data
    if isinstance(data, (bytes, bytearray)): return data.decode('utf-8', errors='ignore')
    return str(data)

def sanitize_notary_text(text: Any) -> str:
    """Institutional Sanitizer: Purges PDF-breaking Unicode."""
    if text is None: return ""
    s = force_str(text)
    
    # Critical replacements for standard PDF fonts
    repls = {
        "\u2122": "(TM)", 
        "\u00ae": "(R)",
        "\u2014": "-", 
        "\u2013": "-", 
        "\u2019": "'", 
        "\u2018": "'",
        "\u201c": '"', 
        "\u201d": '"', 
        "\u2022": "*", 
        "\u2026": "...",
        "α": "alpha", "β": "beta", "γ": "gamma", 
        "Å": "A", "²": "2", "³": "3"
    }
    for char, rep in repls.items():
        s = s.replace(char, rep)
    
    # Nuclear fallback: Strip anything that cannot be encoded in Latin-1
    return s.encode('latin-1', 'ignore').decode('latin-1')
