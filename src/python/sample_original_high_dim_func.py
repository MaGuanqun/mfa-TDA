import numpy as np
import pyvista as pv
import argparse


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


def quartic_potential_2(point):
    """
    Quartic potential function:
    f(x, y, z) = 0.25*(x^4 + y^4) + 0.5*(1 - z)*x^2 + 0.5*cos(z)*y^2
    """
    x, y, z = point
    return 0.25 * (x ** 4 + y ** 4) + 0.5 * (1 - z) * x ** 2 + 0.5 * np.cos(z) * y ** 2

def rotating_quartic_multiwell(point):
    """
    f(x,y,t) = 1/4 * [ (x cos t + y sin t)^2 - 1 ]^2
             + 1/4 * [ (-x sin t + y cos t)^2 - 1 ]^2
    """
    x,y,t= point
    c = np.cos(t)
    s = np.sin(t)
    u =  x * c + y * s
    v = -x * s + y * c
    return 0.25 * ((u**2 - 1)**2 + (v**2 - 1)**2)


parser = argparse.ArgumentParser(description='sample functions.')


parser.add_argument('-i', '--function_name', type=str, default='function_name', help='input function name to sample')
parser.add_argument('-o', '--output_name', type=str, default='file_name.vtk', help='output file to compute critical points tracking')


args = parser.parse_args()


if args.function_name in ['quartic_potential_2', 'rotating_gaussian', 'rotating_quartic_multiwell']:
    min=[-2,-2,0]
    max = [2,2,4]

# Grid definition
nx, ny, nz = 100, 100, 1
x = np.linspace(min[0], max[0], nx)
y = np.linspace(min[1], max[1], ny)
z = np.array([0])

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')  # shape: (nx, ny, nz)

# Flatten the coordinates into a (N, 3) array
points = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)

# Evaluate rotating_gaussian at each point
if args.function_name == 'rotating_gaussian':
    values = np.array([rotating_gaussian(p) for p in points])  # shape (N,)
elif args.function_name == 'quartic_potential_2':
    values = np.array([quartic_potential_2(p) for p in points])
elif args.function_name == 'rotating_quartic_multiwell':
    values = np.array([rotating_quartic_multiwell(p) for p in points])

# Reshape back to grid shape
values = values.reshape((nx, ny, nz))

# Create the pyvista StructuredGrid (Structured = regular topology, but geometry may be warped)
grid = pv.StructuredGrid()
grid.points = points
grid.dimensions = (nx, ny, nz)

# Add scalar field
grid["var0"] = values.ravel(order='F')  # VTK expects Fortran order (z-fastest)

# Save to .vts file (structured grid format)
grid.save(args.output_name)
