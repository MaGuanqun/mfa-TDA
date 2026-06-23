#!/usr/bin/env python3
"""
Extract critical points (zeros) of the discrete gradient vector field stored in a
.vff file, slice by slice in time, and write them as seeds for the Feature Flow
Fields tracker (feature_flow_fields).

Why not TTK/ParaView?
    TTK's critical-point machinery is for *scalar* fields (PL Morse theory); it has
    no native "zero of a vector field" filter. To stay self-consistent with FFF we
    must extract zeros of the *same* multilinearly-interpolated vector field v that
    FFF integrates -- so we do per-cell zero finding directly on the .vff field.
    (A TTK-only approximation would be ScalarFieldCriticalPoints on |v|^2 keeping
    minima with value ~0, but that only snaps to grid vertices and is less
    accurate; see --help notes.)

Method:
    For each requested time slice t_i, find the zeros of the spatial vector field
    v(.,t_i) reconstructed by (bi/tri)linear interpolation inside each grid cell.
      * C == 2 (2D spatial): closed-form bilinear solve (up to 2 roots per cell),
        fully vectorized.
      * C  > 2 (e.g. 3D spatial): sign-change prefilter + Newton per candidate cell.

Output CSV (one row per critical point), with a header line:
    x,y,t,type            (D==3, C==2)
    x,y,z,t,type          (D==4, C==3)
where `type` is a numeric class code (gradient field => symmetric Jacobian):
    0 = minimum (source), 1 = saddle, 2 = maximum (sink), -1 = degenerate.
The FFF tracker only reads the first D columns; `type` is informational.

Usage:
    python3 vff_critical_points.py -i field.vff -o seeds.csv [--t-stride 1] [--tol 1e-9]
"""

import argparse
import sys

import numpy as np

from vff_io import load_vff


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def axis_resolution(vec):
    """vec shape (n[D-1],...,n[0],C) -> n in axis order [n_x, n_y, ..., n_t]."""
    return list(vec.shape[:-1])[::-1]


def spacing_of(dmin, dmax, n_axis):
    n = np.asarray(n_axis, dtype=float)
    return (np.asarray(dmax) - np.asarray(dmin)) / np.maximum(n - 1.0, 1.0)


def classify_2d(jxx, jxy, jyx, jyy):
    """Classify a 2D critical point from its Jacobian (sign-only, robust to scale)."""
    det = jxx * jyy - jxy * jyx
    tr = jxx + jyy
    out = np.full(det.shape, -1, dtype=np.int32)            # degenerate
    out[det < 0] = 1                                        # saddle
    pos = det > 0
    out[pos & (tr > 0)] = 0                                 # minimum / source
    out[pos & (tr < 0)] = 2                                 # maximum / sink
    return out


# ---------------------------------------------------------------------------
# 2D bilinear closed-form extraction (C == 2)
# ---------------------------------------------------------------------------
def extract_slice_bilinear(slc, dmin, spacing, t_value, eps):
    """
    slc      : (ny, nx, 2) vector field on one time slice
    returns list of (x, y, t, type)
    a is the x-parameter in [0,1], b is the y-parameter in [0,1].
    bilinear:  f(a,b) = c0 + c1 a + c2 b + c3 a b
    """
    Fx = slc[:, :, 0]
    Fy = slc[:, :, 1]

    # cell corners: index [j, i] = (y=j, x=i)
    def corners(F):
        return F[:-1, :-1], F[:-1, 1:], F[1:, :-1], F[1:, 1:]   # 00,10,01,11

    fx00, fx10, fx01, fx11 = corners(Fx)
    fy00, fy10, fy01, fy11 = corners(Fy)

    # necessary condition: both components must change sign over the 4 corners
    fx_lo = np.minimum.reduce([fx00, fx10, fx01, fx11])
    fx_hi = np.maximum.reduce([fx00, fx10, fx01, fx11])
    fy_lo = np.minimum.reduce([fy00, fy10, fy01, fy11])
    fy_hi = np.maximum.reduce([fy00, fy10, fy01, fy11])
    cand = (fx_lo <= 0) & (fx_hi >= 0) & (fy_lo <= 0) & (fy_hi >= 0)
    if not cand.any():
        return []

    jj, ii = np.nonzero(cand)            # cell (j=y, i=x) indices

    k0 = fx00[jj, ii]
    k1 = fx10[jj, ii] - k0
    k2 = fx01[jj, ii] - k0
    k3 = fx00[jj, ii] - fx10[jj, ii] - fx01[jj, ii] + fx11[jj, ii]

    l0 = fy00[jj, ii]
    l1 = fy10[jj, ii] - l0
    l2 = fy01[jj, ii] - l0
    l3 = fy00[jj, ii] - fy10[jj, ii] - fy01[jj, ii] + fy11[jj, ii]

    # quadratic in b:  A b^2 + B b + C = 0
    A = l2 * k3 - k2 * l3
    B = (l0 * k3 + l2 * k1) - (k0 * l3 + k2 * l1)
    Cc = l0 * k1 - k0 * l1

    results = []

    def add_roots(b):
        """given candidate b (same shape as ii), recover a and emit valid points."""
        valid = np.isfinite(b) & (b >= -eps) & (b <= 1 + eps)
        if not valid.any():
            return
        bb = np.clip(b[valid], 0.0, 1.0)
        denom = k1[valid] + k3[valid] * bb
        num = -(k0[valid] + k2[valid] * bb)
        ok = np.abs(denom) > eps
        if not ok.any():
            return
        aa = np.full_like(bb, np.nan)
        aa[ok] = num[ok] / denom[ok]
        good = np.isfinite(aa) & (aa >= -eps) & (aa <= 1 + eps)
        if not good.any():
            return
        aa = np.clip(aa[good], 0.0, 1.0)
        bbg = bb[good]
        idx = np.nonzero(valid)[0][good]          # back to candidate-cell index

        ig = ii[idx]
        jg = jj[idx]
        x = dmin[0] + (ig + aa) * spacing[0]
        y = dmin[1] + (jg + bbg) * spacing[1]

        # Jacobian for classification (scaled, sign preserved)
        jxx = (k1[idx] + k3[idx] * bbg) / spacing[0]
        jxy = (k2[idx] + k3[idx] * aa) / spacing[1]
        jyx = (l1[idx] + l3[idx] * bbg) / spacing[0]
        jyy = (l2[idx] + l3[idx] * aa) / spacing[1]
        typ = classify_2d(jxx, jxy, jyx, jyy)

        for xi, yi, ti in zip(x, y, typ):
            results.append((xi, yi, t_value, int(ti)))

    quad = np.abs(A) > eps
    disc = B * B - 4 * A * Cc
    if quad.any():
        with np.errstate(invalid="ignore"):
            sq = np.where(disc >= 0, np.sqrt(np.abs(disc)), np.nan)
        for sign in (+1.0, -1.0):
            b = np.full_like(B, np.nan)
            b[quad] = (-B[quad] + sign * sq[quad]) / (2 * A[quad])
            add_roots(b)

    lin = (~quad) & (np.abs(B) > eps)
    if lin.any():
        b = np.full_like(B, np.nan)
        b[lin] = -Cc[lin] / B[lin]
        add_roots(b)

    return results


# ---------------------------------------------------------------------------
# generic N-d multilinear extraction via Newton (C >= 3)
# ---------------------------------------------------------------------------
def _multilinear_value_jac(u, corner_vals, C):
    """
    u           : (C,) local coords in [0,1]^C
    corner_vals : (2^C, C) vector at each corner; corner k bit d selects u_d vs 1-u_d
    returns f (C,) and J (C, C) = df/du
    """
    ncorner = corner_vals.shape[0]
    f = np.zeros(C)
    J = np.zeros((C, C))
    for k in range(ncorner):
        bits = [(k >> d) & 1 for d in range(C)]
        w = 1.0
        for d in range(C):
            w *= u[d] if bits[d] else (1.0 - u[d])
        f += w * corner_vals[k]
        for d in range(C):
            dw = 1.0
            for e in range(C):
                if e == d:
                    dw *= 1.0 if bits[e] else -1.0
                else:
                    dw *= u[e] if bits[e] else (1.0 - u[e])
            J[:, d] += dw * corner_vals[k]
    return f, J


def extract_slice_newton(slc, dmin, spacing, t_value, eps, C, max_iter=30):
    """
    slc shape (n_{C-1}, ..., n_0, C). Spatial axes only (time already fixed).
    Generic per-cell sign-change prefilter + Newton. dmin/spacing are length C
    (spatial axes only).
    """
    spatial_shape = slc.shape[:-1]                 # (n_{C-1}, ..., n_0)  slow..fast
    # cell grid: each axis has n-1 cells
    cell_counts = [s - 1 for s in spatial_shape]
    if any(c <= 0 for c in cell_counts):
        return []

    # build the 2^C corner offset list in array-index (slow..fast) order
    ncorner = 1 << C
    # axis order in array is reversed vs. "x first"; map local dim d (x=0) to array axis
    # array axis for spatial dim d (x-first) is (C-1-d)
    results = []

    # iterate cells with nested ranges; prefilter via sign change is done per cell
    # (use np.ndindex over cell grid in array order)
    import itertools

    ranges = [range(c) for c in cell_counts]       # array-order (slow..fast)
    for cell in itertools.product(*ranges):
        # gather 2^C corners
        corner_vals = np.empty((ncorner, C))
        all_pos = np.ones(C, dtype=bool)
        all_neg = np.ones(C, dtype=bool)
        for k in range(ncorner):
            # bit b (x-first) selects +1 along x-first dim b
            idx = []
            for ax, base in enumerate(cell):       # ax in array order slow..fast
                d = C - 1 - ax                     # x-first dim
                bit = (k >> d) & 1
                idx.append(base + bit)
            v = slc[tuple(idx)]
            corner_vals[k] = v
            all_pos &= v > 0
            all_neg &= v < 0
        if all_pos.any() or all_neg.any():
            continue                               # some component never changes sign

        # Newton from center
        u = np.full(C, 0.5)
        converged = False
        for _ in range(max_iter):
            f, J = _multilinear_value_jac(u, corner_vals, C)
            try:
                du = np.linalg.solve(J, -f)
            except np.linalg.LinAlgError:
                break
            u = u + du
            if not np.all(np.isfinite(u)):
                break
            u = np.clip(u, -0.25, 1.25)
            if np.linalg.norm(du) < 1e-10:
                converged = True
                break
        if not converged:
            continue
        if np.any(u < -eps) or np.any(u > 1 + eps):
            continue
        u = np.clip(u, 0.0, 1.0)

        # physical position (x-first dims)
        coords = []
        for d in range(C):
            ax = C - 1 - d
            base = cell[ax]
            coords.append(dmin[d] + (base + u[d]) * spacing[d])

        _, J = _multilinear_value_jac(u, corner_vals, C)
        det = np.linalg.det(J)
        typ = 1 if det < 0 else (-1 if abs(det) < eps else 0)
        results.append(tuple(coords) + (t_value, int(typ)))

    return results


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def extract_all(vec, dmin, dmax, t_stride, eps, dedup_tol):
    D = vec.ndim - 1
    C = vec.shape[-1]
    n_axis = axis_resolution(vec)                  # [n_x, n_y, ..., n_t]
    spacing = spacing_of(dmin, dmax, n_axis)
    if C != D - 1:
        print(f"[warn] C ({C}) != D-1 ({D-1}); treating last axis as time, "
              f"first {C} axes as spatial.", file=sys.stderr)

    nt = n_axis[-1]
    t_dmin = dmin[-1]
    t_spacing = spacing[-1]
    spatial_dmin = np.asarray(dmin[:C], dtype=float)
    spatial_spacing = np.asarray(spacing[:C], dtype=float)

    all_rows = []
    for it in range(0, nt, t_stride):
        t_value = t_dmin + it * t_spacing
        slc = vec[it]                              # shape (..., y, x, C); time removed
        if C == 2:
            rows = extract_slice_bilinear(slc, spatial_dmin, spatial_spacing,
                                          t_value, eps)
        else:
            rows = extract_slice_newton(slc, spatial_dmin, spatial_spacing,
                                        t_value, eps, C)
        all_rows.extend(rows)
        print(f"  t-slice {it:4d} (t={t_value:.6g}): {len(rows)} critical points")

    # dedup within rounding tolerance (cells sharing edges may double-report)
    if dedup_tol > 0 and all_rows:
        seen = set()
        uniq = []
        for r in all_rows:
            key = tuple(int(round(c / dedup_tol)) for c in r[:-1])
            if key in seen:
                continue
            seen.add(key)
            uniq.append(r)
        print(f"[dedup] {len(all_rows)} -> {len(uniq)} critical points")
        all_rows = uniq

    return all_rows, D, C


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="input .vff file")
    ap.add_argument("-o", "--output", required=True, help="output seeds .csv file")
    ap.add_argument("--t-stride", type=int, default=1,
                    help="seed from every Nth grid time level "
                         "(default 1 = all grid time levels)")
    ap.add_argument("--tol", type=float, default=1e-9,
                    help="numerical tolerance for solves / membership (default 1e-9)")
    ap.add_argument("--dedup-tol", type=float, default=1e-6,
                    help="merge points closer than this (per axis); 0 disables")
    args = ap.parse_args()

    vec, dmin, dmax = load_vff(args.input)
    D = vec.ndim - 1
    C = vec.shape[-1]
    print(f"[load] {args.input}: D={D} C={C} shape={vec.shape} "
          f"dmin={dmin} dmax={dmax}")

    rows, D, C = extract_all(vec, dmin, dmax, args.t_stride, args.tol, args.dedup_tol)

    # header: spatial axis names + time + type
    axis_names = ["x", "y", "z", "w"][:C]
    header = ",".join(axis_names + ["t", "type"])
    with open(args.output, "w") as f:
        f.write(header + "\n")
        for r in rows:
            f.write(",".join(f"{v:.10g}" for v in r[:-1]) + f",{r[-1]}\n")
    print(f"[write] {args.output}: {len(rows)} critical points "
          f"({len(axis_names)+1} coord cols + type)")


if __name__ == "__main__":
    main()
