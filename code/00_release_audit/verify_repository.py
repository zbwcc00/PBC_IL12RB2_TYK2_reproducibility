from __future__ import annotations

import hashlib
import os
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
MAX_GITHUB_BYTES = 100 * 1024 * 1024
ERRORS = []
IGNORED_PARTS = {".git", "__pycache__", ".mplconfig", ".venv", "work"}


def is_ignored(path: Path) -> bool:
    return any(part in IGNORED_PARTS for part in path.relative_to(REPO).parts)

for path in REPO.rglob("*"):
    if not path.is_file() or is_ignored(path):
        continue
    relative = path.relative_to(REPO).as_posix()
    if not relative.isascii():
        ERRORS.append(f"non-ASCII repository path: {relative}")
    if path.stat().st_size >= MAX_GITHUB_BYTES:
        ERRORS.append(f"file >=100 MiB: {relative}")

text_extensions = {".py", ".r", ".sh", ".ps1", ".md", ".txt", ".yaml", ".yml", ".cff", ".env", ".tsv", ".csv", ".json"}
for path in REPO.rglob("*"):
    if not path.is_file() or is_ignored(path) or path.suffix.lower() not in text_extensions:
        continue
    text = path.read_text(encoding="utf-8", errors="replace")
    if re.search(r"[A-Za-z]:[\\/]Users[\\/]", text, flags=re.I):
        ERRORS.append(f"personal path: {path.relative_to(REPO)}")
    legacy_tokens = ["PBC_" + "PROJECT_ROOT", "PBC_" + "BULK_ROOT", "project" + "_tree"]
    if any(token in text for token in legacy_tokens):
        ERRORS.append(f"legacy path token: {path.relative_to(REPO)}")
    drive_path = re.search(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]", text, flags=re.I)
    mounted_path = re.search("/" + "mnt/" + r"[a-z]/", text, flags=re.I)
    extended_windows_path = ("\\" * 2 + "?" + "\\") in text
    if drive_path or mounted_path or extended_windows_path:
        ERRORS.append(f"absolute local path: {path.relative_to(REPO)}")
    secret_single = re.search(r"(?i)(api[_-]?key|access[_-]?token|secret[_-]?key|password)\s*[:=]\s*'[^']{6,}'", text)
    secret_double = re.search(r'(?i)(api[_-]?key|access[_-]?token|secret[_-]?key|password)\s*[:=]\s*"[^"]{6,}"', text)
    if secret_single or secret_double:
        ERRORS.append(f"potential secret: {path.relative_to(REPO)}")

required = [
    "VERSION",
    "README.md",
    "DATA_AVAILABILITY.md",
    "REPRODUCIBILITY.md",
    "CITATION.cff",
    "LICENSE",
    "manifests/SHA256SUMS.tsv",
    "manifests/figure_asset_map.tsv",
    "manifests/analysis_code_index.tsv",
    "audit/project_code_asset_inventory.tsv",
    "audit/figure_pixel_reproduction_audit.tsv",
    "audit/compact_figure_source_comparison.tsv",
]
for relative in required:
    if not (REPO / relative).exists():
        ERRORS.append(f"missing required file: {relative}")

manifest_path = REPO / "manifests" / "SHA256SUMS.tsv"
if manifest_path.exists():
    manifest = pd.read_csv(manifest_path, sep="	")
    expected_paths = set(manifest["path"].astype(str))
    actual_paths = {
        path.relative_to(REPO).as_posix()
        for path in REPO.rglob("*")
        if path.is_file() and not is_ignored(path) and path != manifest_path
    }
    missing = expected_paths - actual_paths
    unmanifested = actual_paths - expected_paths
    for relative in sorted(missing):
        ERRORS.append(f"manifested file missing: {relative}")
    for relative in sorted(unmanifested):
        ERRORS.append(f"unmanifested file: {relative}")
    for record in manifest.itertuples(index=False):
        path = REPO / record.path
        if not path.exists():
            continue
        digest = hashlib.sha256()
        path_arg = str(path.resolve())
        if os.name == "nt":
            path_arg = chr(92) * 2 + "?" + chr(92) + path_arg
        with open(path_arg, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != record.sha256:
            ERRORS.append(f"SHA256 mismatch: {record.path}")

version_path = REPO / "VERSION"
citation_path = REPO / "CITATION.cff"
if version_path.exists() and citation_path.exists():
    version = version_path.read_text(encoding="utf-8").strip()
    citation = citation_path.read_text(encoding="utf-8")
    if not re.search(rf"(?m)^version:\s*[\"']?{re.escape(version)}[\"']?\s*$", citation):
        ERRORS.append("CITATION.cff version does not match VERSION")

if ERRORS:
    print("REPOSITORY CHECK: FAIL")
    for error in ERRORS:
        print("-", error)
    sys.exit(1)

print("REPOSITORY CHECK: PASS")
print(f"Repository: {REPO}")
print(f"Files: {sum(1 for path in REPO.rglob('*') if path.is_file() and not is_ignored(path))}")
