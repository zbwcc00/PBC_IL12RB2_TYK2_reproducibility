from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
MANIFEST = REPO / "manifests" / "SHA256SUMS.tsv"
IGNORED_PARTS = {".git", "__pycache__", ".mplconfig", ".venv", "work"}


def is_included(path: Path) -> bool:
    if not path.is_file() or path == MANIFEST:
        return False
    return not any(part in IGNORED_PARTS for part in path.relative_to(REPO).parts)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def collect_records() -> list[dict]:
    files = sorted(
        (path for path in REPO.rglob("*") if is_included(path)),
        key=lambda path: path.relative_to(REPO).as_posix(),
    )
    return [
        {
            "path": path.relative_to(REPO).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in files
    ]


records = collect_records()

largest = sorted(records, key=lambda record: record["bytes"], reverse=True)[:20]
pd.DataFrame(largest).to_csv(REPO / "audit" / "largest_files.tsv", sep="\t", index=False, lineterminator="\n")

version = (REPO / "VERSION").read_text(encoding="utf-8").strip()
summary = {
    "release": version,
    "refreshed": date.today().isoformat(),
    "repository": ".",
    "file_count_excluding_checksum_manifest": len(records),
    "total_bytes_excluding_checksum_manifest": sum(record["bytes"] for record in records),
    "largest_file": largest[0]["path"] if largest else None,
    "largest_file_bytes": largest[0]["bytes"] if largest else 0,
    "paths_ascii": all(record["path"].isascii() for record in records),
}
(REPO / "audit" / "release_summary.json").write_text(
    json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
    encoding="utf-8",
)

records = collect_records()
MANIFEST.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(records).to_csv(MANIFEST, sep="\t", index=False, lineterminator="\n")

print(f"Wrote {len(records)} checksum records for release {version}.")
