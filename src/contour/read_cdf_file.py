"""Convert Hurricane Isabel temperature to the raw float32 format for gridded_3d.

From the repository root: python src/contour/read_cdf_file.py
An optional input path and --output path override the defaults.
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
from scipy.ndimage import convolve

# from netCDF4 import Dataset  # Needed only for the commented NetCDF examples.


# file_path = '../../build/examples/boussinesq/boussinesq.nc'
# nc_data = Dataset(file_path, mode='r')


# # print(nc_data)

# u0 = nc_data.variables['u'][1999, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v0 = nc_data.variables['v'][1999, :, :].filled(-2.0)

# u1 = nc_data.variables['u'][2000, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v1 = nc_data.variables['v'][2000, :, :].filled(-2.0)



# max_value = np.max(u1)
# min_value = np.min(u1)

# print(max_value)
# print(min_value)


# velocity_t1500 = np.sqrt(u0**2 + v0**2)
# velocity_t1501 = np.sqrt(u1**2 + v1**2)


# # print(velocity_t1500.shape)

# velocity_t1500.astype('float32').tofile('velocity_t2000.bin')
# velocity_t1501.astype('float32').tofile('velocity_t2001.bin')

# nc_data.close()





# # Replace 'your_file.nc' with the path to your NetCDF file
# file_path = '../../build/examples/vortexstreet/cylinder2d.nc'
# nc_data = Dataset(file_path, mode='r')

# u0 = nc_data.variables['u'][1499, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v0 = nc_data.variables['v'][1499, :, :].filled(-2.0)

# u1 = nc_data.variables['u'][1500, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v1 = nc_data.variables['v'][1500, :, :].filled(-2.0)



# max_value = np.max(u1)
# min_value = np.min(u1)

# print(max_value)
# print(min_value)


# velocity_t1500 = np.sqrt(u0**2 + v0**2)
# velocity_t1501 = np.sqrt(u1**2 + v1**2)


# print(velocity_t1500.shape)

# velocity_t1500.astype('float32').tofile('velocity_t1500.bin')
# velocity_t1501.astype('float32').tofile('velocity_t1501.bin')

# nc_data.close()




# amira_mesh = nib.load('../../build/examples/fluid_simu_ML/3000.am')
# data = amira_mesh.get_fdata()
# print(data.shape)
# print(data)




XDIM = 500
YDIM = 500
ZDIM = 100
DEFAULT_INPUT = (
    Path(__file__).resolve().parent
    / '../../../ori_data/TCf30.bin.gz'
).resolve()
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parent
    / '../../build/examples/hurricane_isabel/TCf30_3d.bin'
).resolve()


def fill_missing_values(data):
    """Fill 1e35/NaN/Inf voxels in place with iterative 3x3x3 valid means.

    Only missing voxels change. Values filled in one pass become neighbors
    for the next pass; points outside the volume do not count as neighbors.
    """
    missing = ~np.isfinite(data) | (data == np.float32(1e35))
    missing_count = int(np.count_nonzero(missing))
    if missing.all():
        raise ValueError('Cannot fill a volume with no valid samples.')

    kernel = np.ones((3, 3, 3), dtype=np.uint8)
    while missing.any():
        neighbor_sum = convolve(
            np.where(missing, 0.0, data), kernel,
            output=np.float64, mode='constant', cval=0.0,
        )
        neighbor_count = convolve(
            (~missing).astype(np.uint8), kernel, mode='constant', cval=0,
        )
        update = missing & (neighbor_count > 0)
        if not update.any():
            raise ValueError('Remaining missing values have no valid neighbors.')
        data[update] = neighbor_sum[update] / neighbor_count[update]
        missing[update] = False

    return missing_count


def convert_hurricane_isabel(input_path, output_path, dims=(XDIM, YDIM, ZDIM)):
    """Read one big-endian gzip volume and write native-endian float32 scalars.

    dims is (x, y, z). C-order shape (z, y, x) preserves the original index
    x + dim_x * (y + dim_y * z), matching Block::read_3d_scalar_data<float>.
    No header, coordinate columns, normalization, or transpose is needed.
    """
    input_path, output_path = Path(input_path), Path(output_path)
    if input_path.resolve() == output_path.resolve():
        raise ValueError('Input and output paths must be different.')
    dim_x, dim_y, dim_z = dims
    if min(dims) <= 0:
        raise ValueError('All grid dimensions must be positive.')

    with gzip.open(input_path, 'rb') as f:
        data = np.frombuffer(f.read(), dtype='>f4')
    expected_size = dim_x * dim_y * dim_z
    if data.size != expected_size:
        raise ValueError(
            f'Data size {data.size} does not match expected size {expected_size}.'
        )

    # astype converts byte order and makes a writable copy; reshape keeps x fastest.
    data = data.astype(np.float32).reshape((dim_z, dim_y, dim_x))
    missing_count = fill_missing_values(data)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data.tofile(output_path)
    print(f'Filled {missing_count} missing voxels using local 3D means.')
    print(f'Wrote {output_path}: {dims} (x, y, z), {data.nbytes} bytes, native float32.')
    print(f'Value range: [{data.min()}, {data.max()}]')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', type=Path, default=DEFAULT_INPUT,
                        help='Input big-endian .bin.gz (default: %(default)s)')
    parser.add_argument('-o', '--output', type=Path, default=DEFAULT_OUTPUT,
                        help='Output native float32 .bin (default: %(default)s)')
    args = parser.parse_args()
    convert_hurricane_isabel(args.input, args.output)




# data = data.reshape((TDIM,ZDIM,YDIM, XDIM))




# slice_z50 = data[0, 50, :, :].copy()

# missing_value = 1e35
# slice_z50[slice_z50 == missing_value] = np.nan

# def nanmean_filter(values):
#     valid_values = values[~np.isnan(values)]
#     if valid_values.size > 0:
#         return np.mean(valid_values)
#     else:
#         return np.nan

# while np.isnan(slice_z50).any():
#     slice_z50 = generic_filter(slice_z50, nanmean_filter, size=3, mode='constant', cval=np.nan)

# slice_z50 = np.nan_to_num(slice_z50, nan=missing_value)

# max_value = np.max(slice_z50)
# min_value = np.min(slice_z50)

# print(max_value)
# print(min_value)


# # print(np.ma.isMaskedArray(slice_z50))

# slice_z50.astype('float32').tofile('../../build/examples/hurricane_isabel/Pf30_50.bin')


# print("Shape of slice_z50:", slice_z50.shape)
# max_value = np.max(slice_z50)
# min_value = np.min(slice_z50)

# print(max_value)
# print(min_value)

# has_nan = np.isnan(slice_z50).any()

# # Print the result
# if has_nan:
#     print("There are NaN values in slice_z50.")
# else:
#     print("There are no NaN values in slice_z50.")


if __name__ == '__main__':
    main()
