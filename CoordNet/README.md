# CoordNet: Data Generation and Visualization Generation for Time-Varying Volumes via a Coordinate-Based Neural Network
Pytorch implementation for CoordNet: Data Generation and Visualization Generation for Time-Varying Volumes via a Coordinate-Based Neural Network

## Prerequisites
- Linux
- CUDA >= 10.0
- Python >= 3.7
- Numpy
- Skimage
- Pytorch >= 1.0

## Data format

The volume at each time step is saved as a .dat file with the little-endian format. The data is stored in column-major order, that is, z-axis goes first, then y-axis, finally x-axis.

## Training models
```
cd Code 
```

- training
```
python3 main.py --mode 'train' --dataset 'Vortex' --applicaion 'temporal'
```

- inference
```
python3 main.py --mode 'inf' --dataset 'Vortex' --application 'temporal'
```

### Fit the vortex street velocity components

From the repository root, use the `siren` Python environment with netCDF4
and CUDA-enabled PyTorch. Run these commands to extract the fields and train
one scalar model per component.

```bash
conda activate siren
python CoordNet/Code/prepare_vortex_street_components.py
bash CoordNet/Code/fit_vortex_street_components.sh
```

The extractor reads `u[1351:1501, :, :]` and `v[1351:1501, :, :]` from
`CoordNet/Data/cylinder2d.nc` (physical time 13.51 through 15.00). It writes
`vortex_street_3d_u.bin` and `vortex_street_3d_v.bin` alongside the source.
Each file contains signed little-endian float32 values in `(time, y, x)`
order with shape `(150, 80, 640)` and x varying fastest.

The trainer runs `u` then `v` using `super-spatial-temporal`, sine activation,
`num_res=5`, 64 initial features, `omega_0=30.0`, and 300 epochs. It uses
the existing vortex street batch size of 16,000 and sampling factor of 1.
Each fit writes `training.log`, its weights (`.pth`), and a float64 TorchScript
model (`-float64.pt`) under `CoordNet/Exp/vortex_street_3d_<component>/`.

CoordNet scales each input axis to `[-1, 1]` and takes coordinates in
`(time, y, x)` order. To convert a normalized model prediction `p` back to
velocity, use `min + (p + 1) * (max - min) / 2` with the min and max of the
corresponding component binary.

## Acknowledgements
This research was supported in part by the start-up fund UDF01002679 of the Chinese University of Hong Kong, Shenzhen, Shenzhen Science and Technology Program ZDSYS20211021111415025, and the U.S. National Science Foundation through grants IIS-1455886, CNS-1629914, DUE-1833129, IIS-1955395, IIS-2101696, and OAC-2104158.

## Citation 
```
@article{han2023coordnet,
  title={CoordNet: Data Generation and Visualization Generation for Time-Varying Volumes via a Coordinate-Based Neural Network},
  author={Han, Jun and Wang, Chaoli},
  journal={IEEE Transactions on Visualization and Computer Graphics},
  volume={29},
  number={12},
  pages={4951--4963},
  year={2023}
}

