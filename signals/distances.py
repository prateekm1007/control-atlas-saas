import numpy as np
from scipy.spatial.distance import cdist

def compute_interface_density(tar_coords, bin_coords, threshold=5.0):
    """Calculates the contact density (rho) as per LAW-153."""
    dists = cdist(tar_coords, bin_coords)
    contacts = np.sum(dists < threshold)
    return contacts

def compute_rmsf(trajectory_coords):
    """Calculates residue-level fluctuation to detect breathing (LAW-150)."""
    mean_pos = np.mean(trajectory_coords, axis=0)
    rmsf = np.sqrt(np.mean(np.sum((trajectory_coords - mean_pos)**2, axis=2), axis=0))
    return rmsf
