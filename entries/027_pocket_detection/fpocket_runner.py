"""
fpocket Subprocess Wrapper
Runs fpocket and parses output files.
"""

import subprocess
import tempfile
import shutil
import os
from pathlib import Path

def run_fpocket(pdb_path: str, timeout: int = 120) -> dict:
    """
    Run fpocket on a PDB file and parse results.
    
    Returns:
    {
        "status": "SUCCESS" | "ERROR",
        "pockets": [
            {
                "pocket_id": str,
                "score": float,
                "druggability": float,
                "volume": float,
                "residues": list,
                "center": [x, y, z]
            },
            ...
        ],
        "output_dir": str | None,
        "error": str | None
    }
    """
    pdb_path = Path(pdb_path).resolve()
    
    if not pdb_path.exists():
        return {
            "status": "ERROR",
            "pockets": [],
            "output_dir": None,
            "error": f"PDB file not found: {pdb_path}"
        }
    
    # Check if fpocket is available
    if shutil.which("fpocket") is None:
        return {
            "status": "ERROR",
            "pockets": [],
            "output_dir": None,
            "error": "fpocket not found in PATH. Install with: conda install -c conda-forge fpocket"
        }
    
    # Create temp working directory
    work_dir = tempfile.mkdtemp(prefix="fpocket_")
    
    try:
        # Copy PDB to work dir (fpocket outputs to same location)
        work_pdb = Path(work_dir) / pdb_path.name
        shutil.copy(pdb_path, work_pdb)
        
        # Run fpocket
        result = subprocess.run(
            ["fpocket", "-f", str(work_pdb)],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=work_dir
        )
        
        if result.returncode != 0:
            return {
                "status": "ERROR",
                "pockets": [],
                "output_dir": None,
                "error": f"fpocket failed: {result.stderr}"
            }
        
        # Parse output
        output_dir = Path(work_dir) / (work_pdb.stem + "_out")
        
        if not output_dir.exists():
            return {
                "status": "ERROR",
                "pockets": [],
                "output_dir": None,
                "error": "fpocket output directory not found"
            }
        
        # Parse pocket info file
        pockets = parse_fpocket_info(output_dir)
        
        return {
            "status": "SUCCESS",
            "pockets": pockets,
            "output_dir": str(output_dir),
            "error": None
        }
        
    except subprocess.TimeoutExpired:
        return {
            "status": "ERROR",
            "pockets": [],
            "output_dir": None,
            "error": f"fpocket timed out after {timeout}s"
        }
    except Exception as e:
        return {
            "status": "ERROR",
            "pockets": [],
            "output_dir": None,
            "error": str(e)
        }


def parse_fpocket_info(output_dir: Path) -> list:
    """Parse fpocket output files to extract pocket data."""
    pockets = []
    output_dir = Path(output_dir)
    
    # Look for pocket PDB files
    pocket_files = sorted(output_dir.glob("pockets/pocket*_atm.pdb"))
    
    # Parse info file for scores
    info_file = output_dir / (output_dir.name.replace("_out", "") + "_info.txt")
    pocket_scores = {}
    
    if info_file.exists():
        pocket_scores = parse_info_file(info_file)
    
    for i, pocket_file in enumerate(pocket_files, 1):
        pocket_id = f"pocket_{i}"
        
        # Extract residues and center from pocket PDB
        residues, center, volume_approx = parse_pocket_pdb(pocket_file)
        
        # Get scores from info file
        scores = pocket_scores.get(i, {})
        
        pockets.append({
            "pocket_id": pocket_id,
            "score": scores.get("score", 0.0),
            "druggability": scores.get("druggability", 0.0),
            "volume": scores.get("volume", volume_approx),
            "residues": residues,
            "center": center,
            "pocket_file": str(pocket_file)
        })
    
    # Sort by score descending
    pockets.sort(key=lambda x: x["score"], reverse=True)
    
    return pockets


def parse_info_file(info_path: Path) -> dict:
    """Parse fpocket info.txt file for pocket scores."""
    scores = {}
    current_pocket = None
    
    try:
        with open(info_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("Pocket"):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            current_pocket = int(parts[1].rstrip(":"))
                            scores[current_pocket] = {}
                        except ValueError:
                            continue
                elif current_pocket and ":" in line:
                    key, _, value = line.partition(":")
                    key = key.strip().lower()
                    value = value.strip()
                    
                    if "score" in key and "drug" not in key:
                        try:
                            scores[current_pocket]["score"] = float(value)
                        except ValueError:
                            pass
                    elif "druggability" in key:
                        try:
                            scores[current_pocket]["druggability"] = float(value)
                        except ValueError:
                            pass
                    elif "volume" in key:
                        try:
                            scores[current_pocket]["volume"] = float(value.split()[0])
                        except (ValueError, IndexError):
                            pass
    except Exception:
        pass
    
    return scores


def parse_pocket_pdb(pocket_file: Path) -> tuple:
    """Extract residues, center, and approximate volume from pocket PDB."""
    residues = set()
    coords = []
    
    try:
        with open(pocket_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        res_name = line[17:20].strip()
                        res_num = int(line[22:26].strip())
                        chain = line[21].strip() or "A"
                        
                        residues.add(f"{chain}:{res_name}{res_num}")
                        
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
    except Exception:
        pass
    
    # Calculate center
    if coords:
        center = [
            sum(c[0] for c in coords) / len(coords),
            sum(c[1] for c in coords) / len(coords),
            sum(c[2] for c in coords) / len(coords)
        ]
    else:
        center = [0.0, 0.0, 0.0]
    
    # Approximate volume from coordinate spread
    if len(coords) > 3:
        x_range = max(c[0] for c in coords) - min(c[0] for c in coords)
        y_range = max(c[1] for c in coords) - min(c[1] for c in coords)
        z_range = max(c[2] for c in coords) - min(c[2] for c in coords)
        volume_approx = x_range * y_range * z_range * 0.5  # Rough ellipsoid
    else:
        volume_approx = 0.0
    
    return list(residues), center, volume_approx
