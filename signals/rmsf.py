import numpy as np

def calculate_rmsf(atom_coords_over_time):
    """Measures structural 'breathing' to identify metastable basins."""
    # coords shape: [frames, atoms, 3]
    avg_coords = np.mean(atom_coords_over_time, axis=0)
    diff = atom_coords_over_time - avg_coords
    rmsf = np.sqrt(np.mean(np.sum(diff**2, axis=2), axis=0))
    return rmsf
