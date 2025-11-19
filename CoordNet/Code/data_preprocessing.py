from netCDF4 import Dataset
import numpy as np
import gzip
from scipy.ndimage import convolve,label
# from morseify_the_data import morseify_3d_slices_auto

from morseify_zero_plateau import (
    morseify_zero_plateau_3d,
    count_changed_entries,
)

# import numpy as np
# import re

# def read_amira_lattice(filename, component=0):
#     with open(filename, 'rb') as f:
#         content = f.read()

#     marker = b"# Data section follows"
#     idx = content.find(marker)
#     if idx < 0:
#         raise ValueError("No data section found")

#     header = content[:idx].decode("ascii", errors="ignore")
#     data_start = idx + len(marker)

#     # Skip newline(s)
#     while content[data_start] in (10, 13):
#         data_start += 1

#     raw = content[data_start:]

#     # Extract dimensions
#     m = re.search(r"define\s+Lattice\s+(\d+)\s+(\d+)\s+(\d+)", header)
#     nx, ny, nz = map(int, m.groups())

#     # Detect dtype
#     if "float" in header:
#         dtype = np.float32
#     elif "double" in header:
#         dtype = np.float64
#     else:
#         raise ValueError("Unsupported datatype")

#     # Detect float[k]
#     comp_match = re.search(r"float\[(\d+)\]", header)
#     n_comp = int(comp_match.group(1)) if comp_match else 1

#     total_vals = nx * ny * nz * n_comp

#     arr = np.frombuffer(raw, dtype=dtype, count=total_vals)
#     arr = arr.reshape((nz, ny, nx, n_comp))

#     return arr[..., component], header, (nx, ny, nz), n_comp


# # -------------------------------------------------------------
# # MAIN CODE: compute magnitude, keep z=1..100, save to binary
# # -------------------------------------------------------------
# filename = "../Data/0000.am"
# output_bin = "../Data/fluid.bin"

# print("Reading U...")
# U, header, (nx, ny, nz), n_comp = read_amira_lattice(filename, component=0)

# print("Reading V...")
# V, _, _, _ = read_amira_lattice(filename, component=1)

# # Compute magnitude
# print("Computing magnitude...")
# mag = np.sqrt(U**2 + V**2).astype(np.float32)

# # -------------------------------------------------------------
# # KEEP ONLY Z = 1..100 (Python index 0..99)
# # -------------------------------------------------------------
# mag_slice = mag[100:200, :, :]   # shape = (100, 512, 512)

# print("Saving z=1..100 to binary file...")
# mag_slice.tofile(output_bin)


def morseify_slice_debug(f2d, mask):
    f = f2d.copy()
    # Just set plateau region to a constant number, e.g., 0.12345
    f[mask] = 0.12345
    return f

save_name = '../Data/cylinder.bin'
file_path = '../Data/pipedcylinder2d.nc'
nc_data = Dataset(file_path, mode='r')


print(nc_data.variables)
print(nc_data.variables['u'].shape)
u0 = nc_data.variables['u'][1400:1500, 78:, 154:].filled(-2.0)  # Shape: (tdim, ydim, xdim)
v0 = nc_data.variables['v'][1400:1500, 78:, 154:].filled(-2.0)
# u0 = nc_data.variables['u'][1400:1500, 78:, 275:].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v0 = nc_data.variables['v'][1400:1500, 78:, 275:].filled(-2.0)
# u1 = nc_data.variables['u'][2000, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v1 = nc_data.variables['v'][2000, :, :].filled(-2.0)



# max_value = np.max(u1)
# min_value = np.min(u1)

# print(max_value)
# print(min_value)


velocity = np.sqrt(u0**2 + v0**2)

# velocity_morse = morseify_zero_plateau_3d(
#     velocity,
#     axis=0,          # slice along time
#     zero_thr=1e-12,  # plateau = |value| < zero_thr
#     frac_eps=0.1,    # strong modification: 10% of range
# )

# changed = count_changed_entries(velocity, velocity_morse)
# print("Changed entries:", changed)

# print("After min/max:", velocity_morse.min(), velocity_morse.max())

print(velocity.shape)
# # velocity_t1501 = np.sqrt(u1**2 + v1**2)

print(np.max(velocity))
print(np.min(velocity))

# max_value = np.max(velocity_morse)
# min_value = np.min(velocity_morse)
# print("max",max_value)
# print("min",min_value)
# print(velocity_morse.shape)



# k = 0
# before = velocity[k]
# plateau_mask = np.abs(before - before.min()) < 1e-6 * (before.max() - before.min())
# after = morseify_slice_debug(before, plateau_mask)

# print("Unique values in plateau before:", np.unique(before[plateau_mask]))
# print("Unique values in plateau after:",  np.unique(after[plateau_mask]))


velocity.astype('float32').tofile(save_name)

# velocity_t1501.astype('float32').tofile('velocity_t2001.bin')

nc_data.close()



# save_name = '../Data/boussinesq_3d.bin'
# file_path = '../Data/boussinesq.nc'
# nc_data = Dataset(file_path, mode='r')

# #(2001, 450, 150)
# print(nc_data.variables)
# print(nc_data.variables['u'].shape)
# u0 = nc_data.variables['u'][300:500, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v0 = nc_data.variables['v'][300:500, :, :].filled(-2.0)

# # u1 = nc_data.variables['u'][2000, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# # v1 = nc_data.variables['v'][2000, :, :].filled(-2.0)



# # max_value = np.max(u1)
# # min_value = np.min(u1)

# # print(max_value)
# # print(min_value)


# velocity = np.sqrt(u0**2 + v0**2)

# print(velocity.shape)
# # # velocity_t1501 = np.sqrt(u1**2 + v1**2)

# max_value = np.max(velocity)
# min_value = np.min(velocity)
# print("max",max_value)
# print("min",min_value)
# # # print(velocity_t1500.shape)

# velocity.astype('float32').tofile(save_name)

# # velocity_t1501.astype('float32').tofile('velocity_t2001.bin')

# nc_data.close()





# # Replace 'your_file.nc' with the path to your NetCDF file
# save_name = '../Data/vortex_street_3d_old.bin' #1351:1501
# file_path = '../Data/cylinder2d.nc'

# nc_data = Dataset(file_path, mode='r')

# print(nc_data.variables)

# u0 = nc_data.variables['u'][1351:1501, :, :].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# v0 = nc_data.variables['v'][1351:1501, :, :].filled(-2.0)

# print(u0.shape)

# # u0 = nc_data.variables['u'][1451:1501, :, 51:151].filled(-2.0)  # Shape: (tdim, ydim, xdim)
# # v0 = nc_data.variables['v'][1451:1501, :, 51:151].filled(-2.0)


# velocity = np.sqrt(u0**2 + v0**2)
# # Normalize the velocity
# velocity_max = np.max(velocity)
# velocity_min = np.min(velocity)
# # velocity_normalized = 2*(velocity - velocity_min) / (velocity_max - velocity_min)-1.0

# # Rearranging the dimensions to ensure last dim changes faster
# # velocity = np.transpose(velocity, (2, 1, 0))  # Change the order of dimensions
# # velocity_max = np.max(velocity_normalized)
# # velocity_min = np.min(velocity_normalized)
# print("max",velocity_max)
# print("min",velocity_min)
# print(velocity.shape)

# velocity.astype('float32').tofile(save_name)
# print("size",velocity.size)

# nc_data.close()

# def fill_nans_local_mean_3d(arr, max_iter=50, kernel_size=3):
#     """
#     Iteratively fill NaNs in a 3D array using the mean of valid neighbors
#     in a local window (default 3x3x3).

#     Parameters
#     ----------
#     arr : np.ndarray
#         3D array with NaNs to fill.
#     max_iter : int
#         Maximum number of iterations to try.
#     kernel_size : int
#         Size of the cubic neighborhood (must be odd).

#     Returns
#     -------
#     filled : np.ndarray
#         New array with NaNs filled where possible.
#     """
#     if arr.ndim != 3:
#         raise ValueError(f"Expected a 3D array, got shape {arr.shape}")

#     if kernel_size % 2 == 0:
#         raise ValueError("kernel_size must be odd (e.g., 3, 5, 7).")

#     # Work on a copy to avoid modifying input in-place
#     filled = arr.copy()

#     # Cubic kernel of ones, e.g. 3x3x3
#     kernel = np.ones((kernel_size, kernel_size, kernel_size), dtype=float)

#     for _ in range(max_iter):
#         nan_mask = np.isnan(filled)
#         if not nan_mask.any():
#             break  # All filled

#         # Replace NaNs with 0 for sum computation
#         filled_zero = np.nan_to_num(filled, nan=0.0)

#         # Sum of neighbors in the local window
#         neighbor_sum = convolve(
#             filled_zero,
#             kernel,
#             mode="constant",
#             cval=0.0
#         )

#         # Number of valid (non-NaN) neighbors
#         valid_mask = (~np.isnan(filled)).astype(float)
#         neighbor_count = convolve(
#             valid_mask,
#             kernel,
#             mode="constant",
#             cval=0.0
#         )

#         # Compute local mean where there is at least one valid neighbor
#         with np.errstate(invalid="ignore", divide="ignore"):
#             local_mean = neighbor_sum / np.maximum(neighbor_count, 1e-12)

#         # We only update NaNs that actually have valid neighbors
#         update_mask = nan_mask & (neighbor_count > 0)

#         if not np.any(update_mask):
#             # Remaining NaNs are completely isolated – cannot be filled
#             break

#         filled[update_mask] = local_mean[update_mask]

#     return filled


# save_name = '../Data/sst.bin'
# file_path = '../Data/sst.mon.mean.nc'
# nc_data = Dataset(file_path, mode='r')

# # --- Load OLR as a 3D array: (time, lat, lon) ---
# # This is usually a masked array; we keep missing values as NaN]
# # print(nc_data.variables['sst'].shape)
# olr_var = nc_data.variables['sst'][1810:2110:,:]     # float32 olr(time, lat, lon)

# olr = olr_var[:].filled(np.nan)            # shape: (tdim, ydim, xdim)
# print("olr shape:", olr.shape)
# print(np.isnan(olr).any())

# # fill_nans_local_mean_3d(olr, max_iter=100, kernel_size=3)
# olr = np.nan_to_num(olr, nan=3)
# # # --- Optional: normalize (ignoring NaNs) ---
# # valid_mask = ~np.isnan(olr)
# # olr_min = np.nanmin(olr)
# # olr_max = np.nanmax(olr)
# # print("olr max:", olr_max)
# # print("olr min:", olr_min)

# # Example: scale to [-1, 1]
# # olr_norm = (2.0 * (olr - olr_min) / (olr_max - olr_min)) - 1.0

# # Put some sentinel for missing values, e.g. -2.0
# # olr_norm[~valid_mask] = -2.0

# # If you don’t need normalization, just do:
# # olr_norm = olr.astype('float32')
# # olr_norm[~valid_mask] = -2.0

# # --- (Optional) rearrange dimensions ---
# # If your downstream code expects last index to be x fastest, you can keep
# # (time, lat, lon) as is. If you wanted (lon, lat, time), you’d do:
# # olr_norm = np.transpose(olr_norm, (2, 1, 0))

# # --- Save to binary ---
# olr.astype('float32').tofile(save_name)
# print("size:", olr.size)

# nc_data.close()



# amira_mesh = nib.load('../../build/examples/fluid_simu_ML/3000.am')
# data = amira_mesh.get_fdata()
# print(data.shape)
# print(data)




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
