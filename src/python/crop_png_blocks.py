#!/usr/bin/env python3
"""
Crop three fixed-size blocks from PNG screenshots for all tracking methods.

Given a dataset and representation (mfa / INR), reads every method figure from
one folder and writes three cropped blocks back into the same folder.

Input filenames (repr comes from the folder, not the filename):
    boussinesq_our_remove_obstacle_high_reso.png
    boussinesq_ttk_remove_obstacle_high_reso.png
    vortex_street_sfff_remove_obstacle_high_reso.png

Output filenames (same folder):
    boussinesq_3d_mfa_our_block1.png
    boussinesq_3d_mfa_our_block2.png
    boussinesq_3d_mfa_our_block3.png

Example:
    python src/python/crop_png_blocks.py \\
        --figs-folder figures/mfa \\
        --dataset boussinesq_3d \\
        --repr mfa
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, List, Sequence, Tuple

from PIL import Image

# ParaView high-res screenshots can exceed PIL's default ~179 MP safety limit.
Image.MAX_IMAGE_PIXELS = None

# ---------------------------------------------------------------------------
# Set block size once here (pixels). All three crops share this size.
# ---------------------------------------------------------------------------
BLOCK_WIDTH = 406 #1300, 1600 (406, vortex_street_3d, INR) (230,vortex_street_3d, mfa) (500, boussinesq_3d, INR)
BLOCK_HEIGHT = 406 #1300, 1600

# Dataset name -> filename prefix used in PNG names.
DATASET_PREFIX = {
    "vortex_street_3d": "vortex_street",
    "boussinesq_3d": "boussinesq",
}

# Suffix shared by all input figures (after "{prefix}_{method}_").
INPUT_NAME_SUFFIX = "remove_obstacle_high_reso.png"

# ---------------------------------------------------------------------------
# Start pixel (left, top) for each of the three blocks.
# Keys: (dataset, repr)  ->  [(x0, y0), (x1, y1), (x2, y2)]
#
# The same three crops are applied to every method (ttk, sfff, our).
# Edit the coordinates below for your screenshots.
# ---------------------------------------------------------------------------
CROP_STARTS: Dict[Tuple[str, str], List[Tuple[int, int]]] = {
    ("vortex_street_3d", "mfa"): [(1625, 2100), (1446, 3225), (1350, 2435)],
    ("vortex_street_3d", "INR"): [(970, 2727), (760, 2274), (1245, 3174)],
    ("boussinesq_3d", "mfa"): [(5514, 3093), (4598, 2260), (3964, 1237)],
    ("boussinesq_3d", "INR"): [(2792, 2039), (5466, 3426), (6587, 1813)],
}

DATASETS = tuple(DATASET_PREFIX.keys())
REPRS = ("mfa", "INR")
METHODS = ("ttk", "sfff", "our")


def normalize_repr(value: str) -> str:
    lowered = value.lower()
    if lowered == "mfa":
        return "mfa"
    if lowered == "inr":
        return "INR"
    raise argparse.ArgumentTypeError(f"repr must be 'mfa' or 'INR', got {value!r}")


def input_png_name(dataset: str, method: str) -> str:
    prefix = DATASET_PREFIX[dataset]
    return f"{prefix}_{method}_{INPUT_NAME_SUFFIX}"


def output_png_name(
    dataset: str,
    repr_name: str,
    method: str,
    block_index: int,
) -> str:
    return f"{dataset}_{repr_name}_{method}_block{block_index}.png"


def get_crop_starts(dataset: str, repr_name: str) -> List[Tuple[int, int]]:
    key = (dataset, repr_name)
    if key not in CROP_STARTS:
        raise KeyError(
            f"No crop configuration for {key}. "
            f"Add an entry to CROP_STARTS in crop_png_blocks.py."
        )
    starts = CROP_STARTS[key]
    if len(starts) != 3:
        raise ValueError(f"Expected 3 start pixels for {key}, got {len(starts)}.")
    return starts


def crop_box(start: Tuple[int, int], width: int, height: int) -> Tuple[int, int, int, int]:
    left, top = start
    return left, top, left + width, top + height


def validate_crop_box(
    box: Tuple[int, int, int, int],
    image_size: Tuple[int, int],
    label: str,
) -> None:
    left, top, right, bottom = box
    img_w, img_h = image_size
    if left < 0 or top < 0:
        raise ValueError(f"{label}: start pixel ({left}, {top}) is negative.")
    if right > img_w or bottom > img_h:
        raise ValueError(
            f"{label}: crop box {box} exceeds image size {image_size} "
            f"(width={img_w}, height={img_h})."
        )
    if right <= left or bottom <= top:
        raise ValueError(f"{label}: invalid crop box {box}.")


def crop_png_blocks(
    input_png: str,
    output_dir: str,
    dataset: str,
    repr_name: str,
    method: str,
    starts: Sequence[Tuple[int, int]],
    block_width: int = BLOCK_WIDTH,
    block_height: int = BLOCK_HEIGHT,
) -> List[str]:
    key = (dataset, repr_name)
    os.makedirs(output_dir, exist_ok=True)

    image = Image.open(input_png)
    saved_paths: List[str] = []

    for index, start in enumerate(starts, start=1):
        box = crop_box(start, block_width, block_height)
        validate_crop_box(box, image.size, f"block {index} of {key}")

        cropped = image.crop(box)
        output_name = output_png_name(dataset, repr_name, method, index)
        output_path = os.path.join(output_dir, output_name)
        cropped.save(output_path)
        saved_paths.append(output_path)
        print(
            f"{method} block {index}: start=({start[0]}, {start[1]}), "
            f"size=({block_width}, {block_height}) -> {output_path}"
        )

    return saved_paths


def crop_all_methods(
    figs_folder: str,
    dataset: str,
    repr_name: str,
    block_width: int = BLOCK_WIDTH,
    block_height: int = BLOCK_HEIGHT,
) -> List[str]:
    figs_folder = os.path.abspath(figs_folder)
    if not os.path.isdir(figs_folder):
        raise FileNotFoundError(f"Figures folder not found: {figs_folder}")

    all_saved: List[str] = []
    missing: List[str] = []
    starts = get_crop_starts(dataset, repr_name)

    for method in METHODS:
        filename = input_png_name(dataset, method)
        input_path = os.path.join(figs_folder, filename)
        if not os.path.isfile(input_path):
            missing.append(filename)
            continue

        print(f"\nProcessing {input_path}")
        saved = crop_png_blocks(
            input_png=input_path,
            output_dir=figs_folder,
            dataset=dataset,
            repr_name=repr_name,
            method=method,
            starts=starts,
            block_width=block_width,
            block_height=block_height,
        )
        all_saved.extend(saved)

    if missing:
        print("\nSkipped missing files:")
        for name in missing:
            print(f"  - {name}")

    if not all_saved:
        raise FileNotFoundError(
            f"No input PNGs found in {figs_folder} for dataset={dataset}. "
            f"Expected names like {input_png_name(dataset, 'our')}."
        )

    print(f"\nSaved {len(all_saved)} cropped image(s) to {figs_folder}")
    return all_saved


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Crop three equal-size blocks from all method PNGs in one folder."
        )
    )
    parser.add_argument(
        "--figs-folder",
        required=True,
        help=(
            "Folder containing input PNGs and where cropped blocks are saved. "
            "Use separate folders for mfa and INR."
        ),
    )
    parser.add_argument(
        "--dataset",
        required=True,
        choices=DATASETS,
        help="Dataset name: vortex_street_3d or boussinesq_3d.",
    )
    parser.add_argument(
        "--repr",
        required=True,
        type=normalize_repr,
        help="Representation: mfa or INR (selects crop coordinates).",
    )
    parser.add_argument(
        "--block-width",
        type=int,
        default=BLOCK_WIDTH,
        help=f"Crop width in pixels (default: {BLOCK_WIDTH}).",
    )
    parser.add_argument(
        "--block-height",
        type=int,
        default=BLOCK_HEIGHT,
        help=f"Crop height in pixels (default: {BLOCK_HEIGHT}).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    crop_all_methods(
        figs_folder=args.figs_folder,
        dataset=args.dataset,
        repr_name=args.repr,
        block_width=args.block_width,
        block_height=args.block_height,
    )


if __name__ == "__main__":
    main()
