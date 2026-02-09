# Deterministic scaffold classification
# This file defines WHAT a scaffold IS, not what it does.

SCAFFOLD_CLASSES = {
    "trp_cage": {
        "geometry": "globular",
        "aspect_ratio": 1.0,
        "radius_A": 12
    },
    "zinc_finger": {
        "geometry": "globular",
        "aspect_ratio": 1.2,
        "radius_A": 14
    },
    "sh3": {
        "geometry": "globular",
        "aspect_ratio": 1.1,
        "radius_A": 15
    },
    "helical": {
        "geometry": "elongated",
        "aspect_ratio": 3.0,
        "radius_A": 4
    },
    "beta_hairpin": {
        "geometry": "elongated",
        "aspect_ratio": 2.5,
        "radius_A": 5
    }
}

def classify_scaffold(name):
    return SCAFFOLD_CLASSES.get(name, None)
