#!/usr/bin/env python3
"""
Convert a matched 4D CSV to 3D by adding 3*t to XYZ.

Input columns expected:
  PositionX, PositionY, PositionZ, t, CriticalType (or CriticalTpe)

Output columns:
  PositionX, PositionY, PositionZ, t, CriticalType (or CriticalTpe)
  RegionId and ColorId are included when present in the input.
"""

import argparse
import csv
import os


def get_critical_column(fieldnames):
    if "CriticalType" in fieldnames:
        return "CriticalType"
    if "CriticalTpe" in fieldnames:
        return "CriticalTpe"
    raise ValueError("Input CSV must contain 'CriticalType' (or 'CriticalTpe').")


def convert_csv(input_csv, output_csv):
    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(input_csv, "r", newline="") as fin:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            raise ValueError("Input CSV has no header.")

        required = ["PositionX", "PositionY", "PositionZ", "t"]
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        crit_col = get_critical_column(reader.fieldnames)
        has_region_id = "RegionId" in reader.fieldnames
        has_color_id = "ColorId" in reader.fieldnames

        out_fields = ["PositionX", "PositionY", "PositionZ", "t", crit_col]
        if has_region_id:
            out_fields.append("RegionId")
        if has_color_id:
            out_fields.append("ColorId")

        with open(output_csv, "w", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fields)
            writer.writeheader()

            for row in reader:
                t_val = float(row["t"])
                shift = 3.0 * t_val
                out_row = {
                    "PositionX": float(row["PositionX"]) + shift,
                    "PositionY": float(row["PositionY"]) + shift,
                    "PositionZ": float(row["PositionZ"]) + shift,
                    "t": t_val,
                    crit_col: row[crit_col],
                }
                if has_region_id:
                    out_row["RegionId"] = row["RegionId"]
                if has_color_id:
                    out_row["ColorId"] = row["ColorId"]
                writer.writerow(out_row)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add 3*t to XYZ; keep t, CriticalType, RegionId, ColorId (if present)."
    )
    parser.add_argument("-i", "--input-csv", required=True, help="Input CSV file")
    parser.add_argument("-o", "--output-csv", required=True, help="Output CSV file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert_csv(args.input_csv, args.output_csv)
    print(f"Wrote: {args.output_csv}")
