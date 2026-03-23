#!/usr/bin/env pvpython
"""
Read a VTI file, tetrahedralize it, extract the isosurface at scalar value 0,
and save the result as a VTP file.
"""

import argparse
import os

from paraview.simple import Contour, Tetrahedralize, XMLImageDataReader, SaveData


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract scalars=0 isosurface from VTI and write VTP."
    )
    parser.add_argument("-i", "--input-vti", required=True, help="Input .vti file")
    parser.add_argument("-o", "--output-vtp", required=True, help="Output .vtp file")
    parser.add_argument(
        "--array-name",
        default="scalars",
        help="Point scalar array name for contouring (default: scalars)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Read VTI
    reader = XMLImageDataReader(registrationName="InputVTI", FileName=[args.input_vti])
    reader.UpdatePipeline()

    # Tetrahedralize
    tetra = Tetrahedralize(registrationName="Tetrahedralize1", Input=reader)
    tetra.UpdatePipeline()

    # Contour at isovalue 0
    contour = Contour(registrationName="Contour1", Input=tetra)
    contour.ContourBy = ["POINTS", args.array_name]
    contour.Isosurfaces = [0.0]
    contour.PointMergeMethod = "Uniform Binning"
    contour.UpdatePipeline()

    out_dir = os.path.dirname(args.output_vtp)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    SaveData(args.output_vtp, proxy=contour)
    print(f"Wrote: {args.output_vtp}")


if __name__ == "__main__":
    main()
