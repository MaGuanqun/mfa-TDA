"""
morseify_field_auto.py

Automatically detect and Morse-ify flat plateau regions WITHOUT knowing
the plateau value. Works on any constant or near-constant region.

Requires:
    numpy
    scipy.ndimage
"""

import numpy as np
from scipy.ndimage import label, maximum_filter, minimum_filter


def detect_plateau_mask_auto(f2d: np.ndarray, delta: float) -> np.ndarray:
    """
    Detect plateau regions automatically using local range test.

    plateau_mask[y,x] = True  if  max(local) - min(local) < delta
    """
    # Local neighborhood (5x5). Change size if needed.
    max_local = maximum_filter(f2d, size=5, mode="nearest")
    min_local = minimum_filter(f2d, size=5, mode="nearest")

    local_range = max_local - min_local
    return local_range < delta


def morseify_slice_2d_auto(
    f2d: np.ndarray,
    delta_frac: float = 1e-4,
    frac_eps: float = 0.05,
    component_min_size: int = 20,
    make_minimum: bool = True,
    strategy: str = "replace",
):
    """
    Strong Morse-ification of 2D slice WITHOUT any known plateau value.

    - Finds all flat/near-flat regions automatically
    - Overwrites each plateau with a smooth paraboloid bump/well
    """

    if f2d.ndim != 2:
        raise ValueError("Expecting 2D array.")

    f = f2d.copy()
    data_min = float(f.min())
    data_max = float(f.max())
    data_range = data_max - data_min

    # Plateau detection threshold
    delta = delta_frac * data_range

    # Auto-detect plateau regions
    mask = detect_plateau_mask_auto(f, delta)

    if not mask.any():
        return f

    # Connected plateau components
    labeled, ncomp = label(mask)

    sign = -1.0 if make_minimum else 1.0
    eps = frac_eps * data_range if data_range > 0 else frac_eps

    for comp_id in range(1, ncomp + 1):
        comp_mask = (labeled == comp_id)
        idx = np.argwhere(comp_mask)

        if idx.shape[0] < component_min_size:
            continue

        # center of component
        c = idx.mean(axis=0)

        # radii
        R = np.max(np.abs(idx - c), axis=0).astype(float)
        R[R == 0] = 1.

        # normalized coords
        u = (idx - c) / R
        r2 = np.sum(u**2, axis=1)

        phi = (1.0 - r2)**2
        phi[r2 > 1.0] = 0.0

        for (p, val) in zip(idx, phi):
            if val <= 0:
                continue
            y, x = p

            if strategy == "replace":
                f[y, x] = data_min + sign * eps * val
            else:
                f[y, x] = f[y, x] + sign * eps * val

    return f


def morseify_3d_slices_auto(
    F: np.ndarray,
    axis: int = 0,
    delta_frac: float = 1e-4,
    frac_eps: float = 0.05,
    component_min_size: int = 20,
    make_minimum: bool = True,
    strategy: str = "replace",
):
    """
    Apply automatic plateau detection and Morseification to each 2D slice
    of a 3D array.

    No plateau_value is required. All flat regions are processed.
    """

    if F.ndim != 3:
        raise ValueError("Expecting 3D array.")

    G = np.moveaxis(F, axis, 0).copy()  # (Nslices, Y, X)

    for k in range(G.shape[0]):
        G[k] = morseify_slice_2d_auto(
            G[k],
            delta_frac=delta_frac,
            frac_eps=frac_eps,
            component_min_size=component_min_size,
            make_minimum=make_minimum,
            strategy=strategy,
        )

    return np.moveaxis(G, 0, axis)


