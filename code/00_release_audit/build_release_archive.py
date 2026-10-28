from __future__ import annotations

import argparse
import hashlib
import zipfile
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
ARCHIVE_ROOT = "PBC_IL12RB2_TYK2_reproducibility"
FIXED_TIME = (2026, 10, 28, 0, 0, 0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

manifest_path = REPO / "manifests" / "SHA256SUMS.tsv"
manifest = pd.read_csv(manifest_path, sep="\t")
relative_paths = sorted(set(manifest["path"].astype(str)) | {"manifests/SHA256SUMS.tsv"})

for relative in relative_paths:
    if not relative.isascii():
        raise RuntimeError(f"Non-ASCII archive path: {relative}")
    if not (REPO / relative).is_file():
        raise FileNotFoundError(relative)

args.output.parent.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(args.output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
    for relative in relative_paths:
        source = REPO / relative
        info = zipfile.ZipInfo(f"{ARCHIVE_ROOT}/{relative}", date_time=FIXED_TIME)
        info.compress_type = zipfile.ZIP_DEFLATED
        info.external_attr = 0o100644 << 16
        archive.writestr(info, source.read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)

digest = sha256(args.output)
checksum_path = args.output.with_suffix(args.output.suffix + ".sha256")
checksum_path.write_text(f"{digest}  {args.output.name}\n", encoding="ascii")
print(f"Archive: {args.output}")
print(f"Files: {len(relative_paths)}")
print(f"SHA256: {digest}")
