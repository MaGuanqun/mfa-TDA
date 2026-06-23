"""
Read/write the .vff (vector flow field) binary format used by feature_flow_fields.

Format (little-endian, self-describing):
    magic   char[4]      "VFF1"
    dtype   uint32       0 = float64, 1 = float32   (applies to the data section)
    D       uint32       number of domain axes; LAST AXIS = time
    C       uint32       vector components per node (= D-1 for a gradient field)
    n[D]    uint32 * D   grid resolution per axis; axis 0 (x) varies fastest
    dmin[D] float64 * D  domain minimum per axis
    dmax[D] float64 * D  domain maximum per axis
    data    dtype * (C * prod(n))

Data layout (AoS, components innermost):
    offset = ( i0 + n0*( i1 + n1*( i2 + ... ) ) ) * C + c
i.e. a node's C components are contiguous, then x (fastest), y, ..., time (slowest).

NumPy convenience:
    Build `vec` with shape (n_time, ..., n_y, n_x, C) -- spatial/time axes in REVERSE
    order (slowest first), x last before C -- and pass it to `save_vff`. Internally we
    only need it C-contiguous in that shape, which already matches the on-disk order.
"""

import numpy as np

_MAGIC = b"VFF1"


def save_vff(filename, vec, dmin, dmax, dtype="float64"):
    """
    Save a vector field.

    vec   : np.ndarray of shape (n[D-1], ..., n[1], n[0], C)
            i.e. axes ordered time(slowest) ... y, x(fastest), then components.
            (This is the natural C-order layout matching the on-disk format.)
    dmin  : sequence of length D (domain min per axis, axis 0 = x ... last = time)
    dmax  : sequence of length D
    dtype : "float64" or "float32"
    """
    vec = np.asarray(vec)
    D = vec.ndim - 1
    C = vec.shape[-1]
    # vec.shape = (n[D-1], ..., n[0], C); recover n in axis order (x first)
    n = list(vec.shape[:-1])[::-1]  # -> [n0(x), n1(y), ..., n[D-1](time)]

    dmin = np.asarray(dmin, dtype=np.float64)
    dmax = np.asarray(dmax, dtype=np.float64)
    assert dmin.size == D and dmax.size == D, "dmin/dmax length must equal D"

    np_dtype = np.float32 if dtype == "float32" else np.float64
    dtype_flag = 1 if dtype == "float32" else 0

    with open(filename, "wb") as f:
        f.write(_MAGIC)
        np.array([dtype_flag, D, C], dtype="<u4").tofile(f)
        np.array(n, dtype="<u4").tofile(f)
        dmin.astype("<f8").tofile(f)
        dmax.astype("<f8").tofile(f)
        np.ascontiguousarray(vec, dtype=np_dtype).tofile(f)

    print(f"[save_vff] wrote {filename}: D={D} C={C} n={n} dtype={dtype}")


def load_vff(filename):
    """
    Load a vector field. Returns (vec, dmin, dmax) where vec has shape
    (n[D-1], ..., n[0], C) (time slowest ... x fastest, components last).
    """
    with open(filename, "rb") as f:
        magic = f.read(4)
        if magic != _MAGIC:
            raise ValueError(f"bad magic {magic!r} in {filename}")
        dtype_flag, D, C = np.fromfile(f, dtype="<u4", count=3)
        D = int(D)
        C = int(C)
        n = np.fromfile(f, dtype="<u4", count=D).astype(int)
        dmin = np.fromfile(f, dtype="<f8", count=D)
        dmax = np.fromfile(f, dtype="<f8", count=D)
        np_dtype = np.float32 if dtype_flag == 1 else np.float64
        count = int(C) * int(np.prod(n))
        data = np.fromfile(f, dtype=np_dtype, count=count)

    shape = list(n[::-1]) + [C]  # (time, ..., y, x, C)
    vec = data.reshape(shape)
    return vec, dmin, dmax
