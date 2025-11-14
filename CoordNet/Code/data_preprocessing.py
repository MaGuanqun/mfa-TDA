from netCDF4 import Dataset
import numpy as np
import gzip
from scipy.ndimage import convolve

# import meshio
# 

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





# Replace 'your_file.nc' with the path to your NetCDF file
file_path = '../Data/cylinder2d.nc'
nc_data = Dataset(file_path, mode='r')

u0 = nc_data.variables['u'][1201:1501, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
v0 = nc_data.variables['v'][1201:1501, :, :].filled(-2.0)

# print(u0.shape)

# u0 = nc_data.variables['u'][1451:1501, :, 51:151].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v0 = nc_data.variables['v'][1451:1501, :, 51:151].filled(-2.0)


velocity = np.sqrt(u0**2 + v0**2)
# Normalize the velocity
velocity_max = np.max(velocity)
velocity_min = np.min(velocity)
# velocity_normalized = 2*(velocity - velocity_min) / (velocity_max - velocity_min)-1.0

# Rearranging the dimensions to ensure last dim changes faster
# velocity = np.transpose(velocity, (2, 1, 0))  # Change the order of dimensions
# velocity_max = np.max(velocity_normalized)
# velocity_min = np.min(velocity_normalized)
print("max",velocity_max)
print("min",velocity_min)
print(velocity.shape)

velocity.astype('float32').tofile('../Data/vortex_street_3d.bin')
print("size",velocity.size)

nc_data.close()




# amira_mesh = nib.load('../../build/examples/fluid_simu_ML/3000.am')
# data = amira_mesh.get_fdata()
# print(data.shape)
# print(data)


def fill_nans_local_mean_3d(arr, max_iter=50, kernel_size=3):
    """
    Iteratively fill NaNs in a 3D array using the mean of valid neighbors
    in a local window (default 3x3x3).

    Parameters
    ----------
    arr : np.ndarray
        3D array with NaNs to fill.
    max_iter : int
        Maximum number of iterations to try.
    kernel_size : int
        Size of the cubic neighborhood (must be odd).

    Returns
    -------
    filled : np.ndarray
        New array with NaNs filled where possible.
    """
    if arr.ndim != 3:
        raise ValueError(f"Expected a 3D array, got shape {arr.shape}")

    if kernel_size % 2 == 0:
        raise ValueError("kernel_size must be odd (e.g., 3, 5, 7).")

    # Work on a copy to avoid modifying input in-place
    filled = arr.copy()

    # Cubic kernel of ones, e.g. 3x3x3
    kernel = np.ones((kernel_size, kernel_size, kernel_size), dtype=float)

    for _ in range(max_iter):
        nan_mask = np.isnan(filled)
        if not nan_mask.any():
            break  # All filled

        # Replace NaNs with 0 for sum computation
        filled_zero = np.nan_to_num(filled, nan=0.0)

        # Sum of neighbors in the local window
        neighbor_sum = convolve(
            filled_zero,
            kernel,
            mode="constant",
            cval=0.0
        )

        # Number of valid (non-NaN) neighbors
        valid_mask = (~np.isnan(filled)).astype(float)
        neighbor_count = convolve(
            valid_mask,
            kernel,
            mode="constant",
            cval=0.0
        )

        # Compute local mean where there is at least one valid neighbor
        with np.errstate(invalid="ignore", divide="ignore"):
            local_mean = neighbor_sum / np.maximum(neighbor_count, 1e-12)

        # We only update NaNs that actually have valid neighbors
        update_mask = nan_mask & (neighbor_count > 0)

        if not np.any(update_mask):
            # Remaining NaNs are completely isolated – cannot be filled
            break

        filled[update_mask] = local_mean[update_mask]

    return filled


# XDIM = 500
# YDIM = 500
# ZDIM = 100
# TDIM = 1  # Single time step for each file

# # Load the binary data
# file_path = '../Data/TCf30.bin.gz'



# with gzip.open(file_path, 'rb') as f:
#     # Step 2: Read the data as Big Endian float32
#     data = np.frombuffer(f.read(), dtype='>f4')

# expected_size = XDIM * YDIM * ZDIM * TDIM
# if data.size != expected_size:
#     raise ValueError(f"Data size {data.size} does not match expected size {expected_size}.")




# data = data.reshape((TDIM,ZDIM,YDIM, XDIM))


# slice = data[0, :, :, :].copy()  

# print(data.shape)

# missing_value = 1e35
# slice[slice == missing_value] = np.nan


# slice_filled = fill_nans_local_mean_3d(slice)

# slice_filled = np.nan_to_num(slice_filled, nan=missing_value)

# max_value = np.max(slice_filled)
# min_value = np.min(slice_filled)

# print(max_value)
# print(min_value)


# # print(np.ma.isMaskedArray(slice))

# slice_filled.astype('float32').tofile('../Data/hurricane_isabel.bin')


# print("Shape of slice:", slice_filled.shape)
# max_value = np.max(slice_filled)
# min_value = np.min(slice_filled)

# print(max_value)
# print(min_value)

# has_nan = np.isnan(slice_filled).any()

# # Print the result
# if has_nan:
#     print("There are NaN values in slice.")
# else:
#     print("There are no NaN values in slice.")
