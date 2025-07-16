import numpy as np
import pyvista as pv


def rotating_gaussian(domain_pt):
    cx = 0.0
    cy = 0.0
    r1 = 0.7
    phi1 = 0.0
    omega1 = 0.5 * np.pi
    sigmax1 = 0.5
    sigmay1 = 0.5
    sigmax2 = 0.5
    sigmay2 = 0.5
    A1 = 1.0
    A2 = 1.0

    t = domain_pt[2]
    x1 = cx + r1 * np.cos(omega1 * t + phi1)
    y1 = cy + r1 * np.sin(omega1 * t + phi1)

    dx1 = domain_pt[0] - x1
    dy1 = domain_pt[1] - y1

    f1 = A1 * np.exp(-((dx1 ** 2) / (2 * sigmax1 ** 2) + (dy1 ** 2) / (2 * sigmay1 ** 2)))
    f2 = 0.0

    if 1.0 < t < 3.0:
        x2 = cx + r1 * np.cos(-omega1 * t + np.pi + phi1)
        y2 = cy + r1 * np.sin(-omega1 * t + np.pi + phi1)

        dx2 = domain_pt[0] - x2
        dy2 = domain_pt[1] - y2

        f2 = A2 * np.exp(-((dx2 ** 2) / (2 * sigmax2 ** 2) + (dy2 ** 2) / (2 * sigmay2 ** 2)))

        f1 = max(f1, f2)

    return f1
 
 
# Grid definition
nx, ny, nz = 64, 64, 100
x = np.linspace(-2.0, 2.0, nx)
y = np.linspace(-2.0, 2.0, ny)
z = np.linspace(0.0, 4.0, nz)

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')  # shape: (nx, ny, nz)

# Flatten the coordinates into a (N, 3) array
points = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)

# Evaluate rotating_gaussian at each point
values = np.array([rotating_gaussian(p) for p in points])  # shape (N,)

# Reshape back to grid shape
values = values.reshape((nx, ny, nz))

# Create the pyvista StructuredGrid (Structured = regular topology, but geometry may be warped)
grid = pv.StructuredGrid()
grid.points = points
grid.dimensions = (nx, ny, nz)

# Add scalar field
grid["var0"] = values.ravel(order='F')  # VTK expects Fortran order (z-fastest)

# Save to .vts file (structured grid format)
grid.save("rotating_gaussian_raw.vtk")
