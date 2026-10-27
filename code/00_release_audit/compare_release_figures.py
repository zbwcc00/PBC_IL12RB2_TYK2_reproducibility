from __future__ import annotations

import argparse
import hashlib
import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


REPO = Path(__file__).resolve().parents[2]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def locate_entry(names: list[str], suffix: str) -> str | None:
    normalized = suffix.replace("\\", "/")
    matches = [name for name in names if name.replace("\\", "/").endswith(normalized)]
    if len(matches) > 1:
        raise RuntimeError(f"Multiple reference entries match {suffix}: {matches}")
    return matches[0] if matches else None


parser = argparse.ArgumentParser()
parser.add_argument("--reference-zip", type=Path, required=True)
args = parser.parse_args()

figure_rows: list[dict] = []
source_rows: list[dict] = []

with zipfile.ZipFile(args.reference_zip) as archive:
    names = archive.namelist()

    for current in sorted((REPO / "results" / "figures").rglob("*.png")):
        relative = current.relative_to(REPO).as_posix()
        entry = locate_entry(names, relative)
        current_bytes = current.read_bytes()
        if entry is None:
            figure_rows.append({"figure": relative, "status": "NO_REFERENCE", "current_sha256": sha256_bytes(current_bytes)})
            continue

        reference_bytes = archive.read(entry)
        current_image = np.asarray(Image.open(io.BytesIO(current_bytes)).convert("RGBA"), dtype=np.int16)
        reference_image = np.asarray(Image.open(io.BytesIO(reference_bytes)).convert("RGBA"), dtype=np.int16)
        same_dimensions = current_image.shape == reference_image.shape
        if same_dimensions:
            difference = np.abs(current_image - reference_image)
            pixel_identical = bool(np.array_equal(current_image, reference_image))
            changed_pixel_fraction = float(np.any(difference != 0, axis=2).mean())
            mean_absolute_channel_difference = float(difference.mean())
            maximum_channel_difference = int(difference.max())
        else:
            pixel_identical = False
            changed_pixel_fraction = np.nan
            mean_absolute_channel_difference = np.nan
            maximum_channel_difference = np.nan

        figure_rows.append(
            {
                "figure": relative,
                "status": "PIXEL_IDENTICAL" if pixel_identical else "EXPECTED_RENDER_DIFFERENCE",
                "reference_dimensions": f"{reference_image.shape[1]}x{reference_image.shape[0]}",
                "reproduced_dimensions": f"{current_image.shape[1]}x{current_image.shape[0]}",
                "same_dimensions": same_dimensions,
                "changed_pixel_fraction": changed_pixel_fraction,
                "mean_absolute_channel_difference": mean_absolute_channel_difference,
                "maximum_channel_difference": maximum_channel_difference,
                "reference_sha256": sha256_bytes(reference_bytes),
                "current_sha256": sha256_bytes(current_bytes),
            }
        )

    for current in sorted((REPO / "data" / "derived" / "figures").rglob("*")):
        if not current.is_file():
            continue
        relative = current.relative_to(REPO).as_posix()
        entry = locate_entry(names, relative)
        current_bytes = current.read_bytes()
        reference_bytes = archive.read(entry) if entry else None
        source_rows.append(
            {
                "path": relative,
                "status": "BYTE_IDENTICAL" if reference_bytes == current_bytes else ("NO_REFERENCE" if entry is None else "CHANGED"),
                "reference_sha256": sha256_bytes(reference_bytes) if reference_bytes is not None else "",
                "current_sha256": sha256_bytes(current_bytes),
            }
        )

pd.DataFrame(figure_rows).to_csv(
    REPO / "audit" / "figure_pixel_reproduction_audit.tsv",
    sep="\t",
    index=False,
    lineterminator="\n",
)
pd.DataFrame(source_rows).to_csv(
    REPO / "audit" / "compact_figure_source_comparison.tsv",
    sep="\t",
    index=False,
    lineterminator="\n",
)

print(pd.Series(record["status"] for record in figure_rows).value_counts().to_string())
print(pd.Series(record["status"] for record in source_rows).value_counts().to_string())
