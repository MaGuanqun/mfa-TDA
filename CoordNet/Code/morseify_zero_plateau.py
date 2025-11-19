"""
morseify_zero_plateau.py

Strong Morse-ification for 2D/3D scalar fields that contain
zero-valued plateau regions.

This version:
 - Automatically detects near-zero plateau (|f| < zero_thr)
 - Works on each 2D slice independently
 - Replaces each plateau region with a smooth negative paraboloid well
 - Guarantees a single non-degenerate Morse minimum in each plateau
 - No need to know the plateau value (assumes plateau ≈ 0)
 
Requires:
    numpy
    scipy.ndimage (for label)
"""

import numpy as np
from scipy.ndimage import label


# ================================================================
# 1. Strong 2D Morse-ification: zero plateau → negative well
# ================================================================
def morseify_zero_plateau_slice(
    f2d,
    zero_thr=1e-12,
    frac_eps=0.1,
):
    """
    Strongly modify all cells with |f| < zero_thr in a 2D slice.

    For each connected plateau component:
       - Compute component center
       - Compute radii
       - Build smooth radial function phi = (1 - r^2)^2
       - Overwrite plateau with:  f_new = -eps * phi
         (eps = frac_eps * data_range)

    Parameters
    ----------
    f2d : 2D numpy array (H, W)
    zero_thr : float
        Plateau detection threshold.
        Use larger (1e-8 or 1e-6) if plateau isn't exactly zero.
    frac_eps : float
        Well depth as fraction of data range.
        0.05–0.1 is typical.

    Returns
    -------
    f_out : 2D numpy array
    """
    if f2d.ndim != 2:
        raise ValueError("morseify_zero_plateau_slice: input must be 2D")

    f = f2d.copy()
    data_min = float(f.min())
    data_max = float(f.max())
    data_range = max(data_max - data_min, 1.0)

    # ------------------------------------------------------------
    # Plateau region: values very close to 0
    # ------------------------------------------------------------
    mask = np.abs(f) < zero_thr
    if not mask.any():
        return f  # nothing to modify

    # Connected plateau components
    labeled, ncomp = label(mask)

    # Depth of well
    eps = frac_eps * data_range   # negative
    # print("[DEBUG] data_range =", data_range, "eps =", eps)

    # ------------------------------------------------------------
    # Process each plateau component
    # ------------------------------------------------------------
    for comp_id in range(1, ncomp + 1):
        comp_mask = (labeled == comp_id)
        coords = np.argwhere(comp_mask)   # shape (N, 2), rows (y, x)

        if coords.shape[0] == 0:
            continue

        # Center in index space
        c = coords.mean(axis=0)   # float (cy, cx)

        # Radii (extent)
        R = np.max(np.abs(coords - c), axis=0).astype(float)
        R[R == 0.0] = 1.0    # avoid divide-by-zero

        # Normalized coords
        u = (coords - c) / R     # u_ij ~ [-1, 1]
        r2 = np.sum(u**2, axis=1)

        # Smooth bump/well shape
        phi = (1.0 - r2)**2
        phi[r2 > 1.0] = 0.0

        # Overwrite plateau with negative well centered at c
        for (y, x), val in zip(coords, phi):
            if val > 0.0:
                f[y, x] = -eps * val     # < 0

    return f


# ================================================================
# 2. Apply to 3D field slice-by-slice along chosen axis
# ================================================================
def morseify_zero_plateau_3d(
    F,
    axis=0,
    zero_thr=1e-12,
    frac_eps=0.1,
):
    """
    Apply morseify_zero_plateau_slice to each 2D slice of a 3D array.

    Parameters
    ----------
    F : 3D numpy array
        E.g. shape (T, H, W) with axis=0 for time slicing.
    axis : int
        Slice axis (0, 1, or 2).
    zero_thr : float
        Threshold for |f| < zero_thr.
    frac_eps : float
        Well depth fraction.

    Returns
    -------
    F_out : 3D numpy array
    """
    if F.ndim != 3:
        raise ValueError("morseify_zero_plateau_3d: input must be 3D")

    # Move slice axis to front
    G = np.moveaxis(F, axis, 0).copy()   # G shape: (Nslices, H, W)

    # Process each 2D slice
    for k in range(G.shape[0]):
        G[k] = morseify_zero_plateau_slice(G[k], zero_thr=zero_thr, frac_eps=frac_eps)

    # Move axis back
    return np.moveaxis(G, 0, axis)


# ================================================================
# 3. OPTIONAL DEBUG HELPERS
# ================================================================
def count_changed_entries(A, B):
    """
    Count how many entries differ between A and B.
    """
    return np.count_nonzero(A != B)


def print_local_stats(F, mask):
    """
    Print min/max inside mask region.
    """
    vals = F[mask]
    if vals.size > 0:
        print("Masked region min/max:", vals.min(), vals.max())
    else:
        print("Mask is empty.")
