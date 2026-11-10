"""Check VFF layout and derivatives against scalar finite differences."""
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


def potential(x, y, t):
    return x**4/4 + (1-t)*x*x/2 + y**4/4 + math.cos(t)*y*y/2


with tempfile.TemporaryDirectory(prefix='quartic_gradient_') as work:
    output = Path(work) / 'sample.vff'
    executable = sys.argv[1]
    # Unequal dimensions catch transposed axes and component/layout mistakes.
    subprocess.run([executable, '-x', '7', '-y', '6', '-t', '5', '-o', str(output)], check=True)
    with output.open('rb') as stream:
        assert stream.read(4) == b'VFF1'
        assert struct.unpack('<III', stream.read(12)) == (0, 3, 3)
        assert struct.unpack('<3I', stream.read(12)) == (7, 6, 5)
        assert struct.unpack('<3d', stream.read(24)) == (-2., -2., 0.)
        assert struct.unpack('<3d', stream.read(24)) == (2., 2., 4.)
        values = struct.unpack('<630d', stream.read(630*8))
        assert stream.read() == b''
    for it in range(5):
        for iy in range(6):
            for ix in range(7):
                point = [-2 + ix*4/6, -2 + iy*4/5, float(it)]
                for component in range(3):
                    left, right = point.copy(), point.copy()
                    left[component] -= 1e-5
                    right[component] += 1e-5
                    expected = (potential(*right)-potential(*left))/2e-5
                    actual = values[((it*6+iy)*7+ix)*3+component]
                    assert math.isclose(actual, expected, rel_tol=1e-8, abs_tol=1e-8), (point, component, actual, expected)
    rejected = Path(work) / 'invalid.vff'
    result = subprocess.run([executable, '-x', '1', '-o', str(rejected)], capture_output=True)
    assert result.returncode != 0 and not rejected.exists()
    result = subprocess.run([executable, '-o', str(Path(work)/'missing'/'sample.vff')], capture_output=True)
    assert result.returncode != 0
    print('VFF header, dimensions, x-fastest layout, all derivatives and invalid inputs passed.')
