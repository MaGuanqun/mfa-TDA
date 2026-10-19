"""Run with: python -m unittest discover -s src/contour -p test_read_cdf_file.py"""

import contextlib
import gzip
import io
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

from read_cdf_file import convert_hurricane_isabel, fill_missing_values


class HurricaneConversionTests(unittest.TestCase):
    def test_native_float32_x_fastest_without_header(self):
        # Unequal axis lengths expose transposes; coordinate-coded values expose
        # axis swaps even when the output byte count is correct.
        dims = (4, 3, 2)
        expected = [100 * z + 10 * y + x
                    for z in range(2) for y in range(3) for x in range(4)]
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'pressure.bin.gz'
            output = Path(tmp) / 'converted' / 'pressure.bin'
            with gzip.open(source, 'wb') as stream:
                stream.write(struct.pack('>24f', *expected))
            with contextlib.redirect_stdout(io.StringIO()):
                convert_hurricane_isabel(source, output, dims)
            self.assertEqual(output.read_bytes(), struct.pack('=24f', *expected))

    def test_missing_center_uses_3d_neighbors_and_preserves_valid_samples(self):
        data = np.arange(27, dtype=np.float32).reshape(3, 3, 3)
        data[1, 1, 1] = np.float32(1e35)
        valid = data != np.float32(1e35)
        original = data.copy()
        self.assertEqual(fill_missing_values(data), 1)
        self.assertEqual(data[1, 1, 1], original[valid].mean())
        np.testing.assert_array_equal(data[valid], original[valid])

    def test_missing_corner_does_not_count_outside_volume(self):
        data = np.arange(8, dtype=np.float32).reshape(2, 2, 2)
        data[0, 0, 0] = np.nan
        fill_missing_values(data)
        self.assertEqual(data[0, 0, 0], 4.0)

    def test_multiple_passes_fill_sentinel_nan_and_inf(self):
        data = np.full((5, 5, 5), np.float32(1e35))
        data[0, 0, 0] = 7.0
        data[2, 2, 2] = np.nan
        data[4, 4, 4] = np.inf
        self.assertEqual(fill_missing_values(data), 124)
        np.testing.assert_array_equal(data, 7.0)

    def test_all_missing_fails_instead_of_looping(self):
        with self.assertRaisesRegex(ValueError, 'no valid samples'):
            fill_missing_values(np.full((2, 2, 2), np.float32(1e35)))

    def test_wrong_sample_count_does_not_write_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'truncated.bin.gz'
            output = Path(tmp) / 'pressure.bin'
            with gzip.open(source, 'wb') as stream:
                stream.write(struct.pack('>f', 1.0))
            with self.assertRaisesRegex(ValueError, 'expected size 24'):
                convert_hurricane_isabel(source, output, (4, 3, 2))
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
