"""
Toscanini Phase B2 — Protocol Selector
Determines optimal refinement protocol based on failing laws.

Rules:
- LAW-130 (Clashes) + LAW-125 (Rama) + LAW-150 (Rotamer) → Rosetta FastRelax
- LAW-182 (Hydrophobic burial) + LAW-200 (Packing) → OpenMM equilibration
- LAW-110 (Backbone gaps) → Rosetta KIC loop modeling
- Multiple violations → Both (Rosetta first, then OpenMM)
- Default → OpenMM (safer, no license required)
"""
from typing import List, Dict

# Laws best addressed by Rosetta
ROSETTA_LAWS = {"LAW-100", "LAW-120", "LAW-125", "LAW-130", "LAW-135",
                "LAW-145", "LAW-150", "LAW-195"}

# Laws best addressed by OpenMM MD
OPENMM_LAWS  = {"LAW-182", "LAW-200", "LAW-155"}

# Laws requiring loop modeling
LOOP_LAWS    = {"LAW-110", "LAW-160"}

# Protocol options
PROTOCOL_ROSETTA = "rosetta"
PROTOCOL_OPENMM  = "openmm"
PROTOCOL_BOTH    = "both"
PROTOCOL_LOOP    = "loop"


def select_protocol(failing_laws: List[str]) -> Dict:
    """
    Select optimal refinement protocol for given failing laws.

    Args:
        failing_laws: List of failing law IDs (e.g. ["LAW-125", "LAW-150"])

    Returns:
        dict with protocol, rationale, estimated_minutes
    """
    if not failing_laws:
        return {
            "protocol":          PROTOCOL_OPENMM,
            "rationale":         "No violations — light equilibration recommended",
            "estimated_minutes": 3
        }

    failing_set   = set(failing_laws)
    needs_rosetta = bool(failing_set & ROSETTA_LAWS)
    needs_openmm  = bool(failing_set & OPENMM_LAWS)
    needs_loop    = bool(failing_set & LOOP_LAWS)

    # Loop modeling takes priority if backbone gaps
    if needs_loop and needs_rosetta:
        return {
            "protocol":          PROTOCOL_BOTH,
            "rationale":         "Backbone gaps + geometry violations: loop modeling then relax",
            "estimated_minutes": 20,
            "order":             ["loop", "rosetta"]
        }

    if needs_loop:
        return {
            "protocol":          PROTOCOL_LOOP,
            "rationale":         "Backbone gaps detected — KIC loop modeling required",
            "estimated_minutes": 15
        }

    if needs_rosetta and needs_openmm:
        return {
            "protocol":          PROTOCOL_BOTH,
            "rationale":         "Geometry + packing violations: FastRelax then MD equilibration",
            "estimated_minutes": 12,
            "order":             ["rosetta", "openmm"]
        }

    if needs_rosetta:
        return {
            "protocol":          PROTOCOL_ROSETTA,
            "rationale":         "Geometry violations (rama/rotamer/clash): Rosetta FastRelax",
            "estimated_minutes": 5
        }

    if needs_openmm:
        return {
            "protocol":          PROTOCOL_OPENMM,
            "rationale":         "Packing/burial violations: OpenMM equilibration",
            "estimated_minutes": 8
        }

    # Fallback
    return {
        "protocol":          PROTOCOL_OPENMM,
        "rationale":         "Unknown violation type — default OpenMM equilibration",
        "estimated_minutes": 5
    }


def get_openmm_config(failing_laws: List[str]) -> Dict:
    """
    Generate OpenMM protocol config tuned for specific violations.

    Args:
        failing_laws: List of failing law IDs

    Returns:
        Protocol config dict for execute_openmm task
    """
    failing_set = set(failing_laws)

    # More steps for packing/burial issues
    if failing_set & {"LAW-182", "LAW-200"}:
        sim_steps = 2500000  # 5ns
    else:
        sim_steps = 1000000  # 2ns default

    return {
        "temperature_k": 300,
        "sim_steps":     sim_steps,
        "timestep_fs":   2,
        "rationale":     f"Tuned for: {', '.join(sorted(failing_set))}"
    }


def get_rosetta_scoreterms(failing_laws: List[str]) -> Dict:
    """
    Generate Rosetta scoreterm weights tuned for specific violations.

    Args:
        failing_laws: List of failing law IDs

    Returns:
        Weight dict for Rosetta XML generation
    """
    failing_set = set(failing_laws)

    return {
        "rama_weight": "2.0" if "LAW-125" in failing_set else "1.0",
        "dun_weight":  "2.0" if "LAW-150" in failing_set else "1.0",
        "rep_weight":  "1.5" if "LAW-130" in failing_set else "0.55",
        "omega_weight":"2.0" if "LAW-135" in failing_set else "1.0",
        "geom_weight": "2.0" if failing_set & {"LAW-100","LAW-120"} else "1.0",
        "dslf_weight": "2.0" if "LAW-195" in failing_set else "1.0",
    }
