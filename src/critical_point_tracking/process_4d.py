import argparse
import csv
import os
from typing import List

import numpy as np


def read_scalar_array(bin_path: str) -> np.ndarray:
    """
    Read a binary file of doubles into a 1D numpy array.
    """
    return np.fromfile(bin_path, dtype=np.float64)


def process_single_csv(
    csv_path: str,
    scalar_value: float,
    output_path: str,
    keep_all_columns: bool = False,
) -> None:
    """
    Read one CSV, add / rename columns as requested, and write to output_path.

    - Rename PositionX -> v0, PositionY -> v1, PositionZ -> v2
    - Add a new column v3 with constant value = scalar_value for every row
    """
    with open(csv_path, "r", newline="") as f_in:
        reader = csv.reader(f_in)
        rows: List[List[str]] = list(reader)

    if not rows:
        # Empty file, just write header with v0..v3 and return
        with open(output_path, "w", newline="") as f_out:
            writer = csv.writer(f_out)
            writer.writerow(["v0", "v1", "v2", "v3"])
        return

    header = rows[0]
    data_rows = rows[1:]

    # Default behavior: drop all columns except v0/v1/v2 (positions) + v3.
    # (This keeps the output compact and avoids dragging extra TTK columns forward.)
    if not keep_all_columns:
        # Support both the original naming (PositionX/Y/Z) and already-renamed (v0/v1/v2).
        def col_index(pos_name: str, fallback_name: str) -> int:
            if pos_name in header:
                return header.index(pos_name)
            if fallback_name in header:
                return header.index(fallback_name)
            raise RuntimeError(f"Cannot find column '{pos_name}' (or fallback '{fallback_name}') in {csv_path}")

        idx_x = col_index("PositionX", "v0")
        idx_y = col_index("PositionY", "v1")
        idx_z = col_index("PositionZ", "v2")

        with open(output_path, "w", newline="") as f_out:
            writer = csv.writer(f_out)
            writer.writerow(["v0", "v1", "v2", "v3"])
            for row in data_rows:
                # Use exactly the position fields; append v3 scalar value.
                writer.writerow([row[idx_x], row[idx_y], row[idx_z], f"{scalar_value}"])
        return

    # keep_all_columns=True: map old names to new names; keep everything else unchanged.
    renamed_header: List[str] = []
    for name in header:
        if name == "PositionX":
            renamed_header.append("v0")
        elif name == "PositionY":
            renamed_header.append("v1")
        elif name == "PositionZ":
            renamed_header.append("v2")
        else:
            renamed_header.append(name)

    # Ensure v3 appears as a new column at the end
    renamed_header.append("v3")

    with open(output_path, "w", newline="") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(renamed_header)

        for row in data_rows:
            # Pad or trim row so we can safely append v3
            # (we keep all existing columns, just add v3 at the end)
            new_row = list(row)
            new_row.append(f"{scalar_value}")
            writer.writerow(new_row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Attach a scalar from a binary result file to per-index CSVs.\n"
            "For index i, reads input pattern with {index} replaced by i, "
            "renames PositionX/Y/Z to v0/v1/v2, and appends v3 with the scalar value."
        )
    )

    parser.add_argument(
        "-b",
        "--bin",
        dest="bin_path",
        required=True,
        help="Path to result.bin (array of double)",
    )
    parser.add_argument(
        "-p",
        "--pattern",
        dest="input_pattern",
        required=True,
        help=(
            "Pattern for input CSVs', "
            "e.g. 'build/src/vortex/vortex.mfa_8_'"
        ),
    )
    parser.add_argument(
        "-o",
        "--output_pattern",
        dest="output_pattern",
        default=None,
        help=(
            "Pattern for output CSVs "
            "Default: same as input_pattern but with '_update' before extension."
        ),
    )
    parser.add_argument(
        "--keep_all_columns",
        action="store_true",
        help="If set, keep all original columns (after renaming PositionX/Y/Z -> v0/v1/v2) and append v3. "
        "Default: only output v0,v1,v2,v3.",
    )

    args = parser.parse_args()

    scalars = read_scalar_array(args.bin_path)
    n = scalars.shape[0]
    
    print("read steps ",len(scalars))

    # Derive default output pattern if not given
    if args.output_pattern is None:
        # Insert '_update' before file extension
        # Example: input_x.csv -> input_x_update.csv

        output_pattern = args.input_pattern+'_update'
    else:
        output_pattern = args.output_pattern

    for i in range(n):
        csv_in = args.input_pattern+str(i)+'.csv'
        csv_out = output_pattern+str(i)+'.csv'

        if not os.path.exists(csv_in):
            # Skip missing CSVs but continue processing others
            print(f"Warning: input CSV not found for index {i}: {csv_in}")
            continue

        scalar_value = float(scalars[i])
        print(f"Processing index {i}: {csv_in} -> {csv_out} (v3={scalar_value})")
        process_single_csv(
            csv_in,
            scalar_value,
            csv_out,
            keep_all_columns=args.keep_all_columns,
        )


if __name__ == "__main__":
    main()