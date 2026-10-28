from __future__ import annotations

import argparse
import math
import urllib.request
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats


REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "single_cell" / "andrews"
OUT = REPO / "work" / "single_cell" / "andrews"
GENE_ID = "ENSG00000081985"
GENE_SYMBOL = "IL12RB2"

ASSETS = {
    "snRNA": {
        "filename": "34702b40-d0ac-4ec2-8871-776ca012efd0.h5ad",
        "url": "https://datasets.cellxgene.cziscience.com/34702b40-d0ac-4ec2-8871-776ca012efd0.h5ad",
        "bytes": 1543810661,
    },
    "scRNA": {
        "filename": "ff713455-e6b7-45b3-b2dd-1c9dd26251f6.h5ad",
        "url": "https://datasets.cellxgene.cziscience.com/ff713455-e6b7-45b3-b2dd-1c9dd26251f6.h5ad",
        "bytes": 1118963144,
    },
}

TARGETS = {
    "scRNA": {
        "lrNK": {"lrNK"},
        "cNK": {"cNK"},
        "NK_all": {"lrNK", "cNK", "CD8T-cNK", "CD3T-lrNK"},
        "T_NK_union": {"lrNK", "cNK", "CD8T-cNK", "CD3T-lrNK", "CD4T", "CD8T", "Tcell"},
        "CD4T": {"CD4T"},
        "CD8T": {"CD8T"},
    },
    "snRNA": {
        "lrNK": {"lrNK"},
        "T_NK_union": {"lrNK", "CD4T", "Tcell-Doublet"},
        "CD4T": {"CD4T"},
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce Andrews et al. external PBC liver donor-level IL12RB2 pseudobulk analyses."
    )
    parser.add_argument("--scrna-h5ad", type=Path, default=RAW / ASSETS["scRNA"]["filename"])
    parser.add_argument("--snrna-h5ad", type=Path, default=RAW / ASSETS["snRNA"]["filename"])
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download the two public CELLxGENE H5AD assets when they are absent.",
    )
    return parser.parse_args()


def download_if_requested(path: Path, modality: str, enabled: bool) -> None:
    if path.exists():
        return
    if not enabled:
        raise FileNotFoundError(
            f"Missing {path}. Run with --download or place the CELLxGENE asset at this relative path."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    urllib.request.urlretrieve(str(ASSETS[modality]["url"]), temporary)
    observed_bytes = temporary.stat().st_size
    expected_bytes = int(ASSETS[modality]["bytes"])
    if observed_bytes != expected_bytes:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"Size check failed for {modality}: expected {expected_bytes}, observed {observed_bytes}"
        )
    temporary.replace(path)


def disease_group(value: object) -> str | None:
    label = str(value).strip().lower()
    if label == "primary biliary cholangitis":
        return "PBC"
    if label == "normal":
        return "Normal"
    return None


def matrix_sum(matrix: object) -> float:
    if sparse.issparse(matrix):
        return float(matrix.sum())
    return float(np.asarray(matrix).sum())


def library_sum_by_cell(matrix: object) -> float:
    per_cell = np.asarray(matrix.sum(axis=1), dtype=np.float64).reshape(-1)
    return float(per_cell.sum())


def aggregate_modality(path: Path, modality: str) -> pd.DataFrame:
    adata = ad.read_h5ad(path, backed="r")
    try:
        if adata.raw is None:
            raise RuntimeError(f"{path.name} does not contain raw counts")
        raw = adata.raw
        raw_var = raw.var
        if GENE_ID in raw_var.index:
            gene_index = int(raw_var.index.get_loc(GENE_ID))
        else:
            matches = np.flatnonzero(raw_var["feature_name"].astype(str).eq(GENE_SYMBOL).to_numpy())
            if len(matches) != 1:
                raise RuntimeError(f"Expected one {GENE_SYMBOL} feature in {path.name}; found {len(matches)}")
            gene_index = int(matches[0])

        obs = adata.obs.copy()
        obs["group"] = obs["disease"].map(disease_group)
        records: list[dict[str, float | int | str]] = []
        for target, labels in TARGETS[modality].items():
            selected = obs["group"].notna() & obs["author_cell_type"].astype(str).isin(labels)
            target_obs = obs.loc[selected, ["group", "donor_id"]]
            for (group, donor_id), donor_cells in target_obs.groupby(["group", "donor_id"], observed=True):
                indices = obs.index.get_indexer(donor_cells.index)
                indices = np.sort(indices[indices >= 0])
                donor_matrix = raw.X[indices, :]
                gene_sum = matrix_sum(donor_matrix[:, gene_index])
                library_sum = library_sum_by_cell(donor_matrix)
                cpm = gene_sum / library_sum * 1e6 if library_sum > 0 else float("nan")
                records.append(
                    {
                        "group": group,
                        "donor_id": donor_id,
                        "n_cells": int(len(indices)),
                        "gene_sum": gene_sum,
                        "lib_sum": library_sum,
                        "CPM": cpm,
                        "logCPM": math.log2(cpm + 0.5),
                        "modality": modality,
                        "target": target,
                    }
                )
        frame = pd.DataFrame(records)
        frame["group_order"] = frame["group"].map({"PBC": 0, "Normal": 1})
        frame["target_order"] = frame["target"].map({name: index for index, name in enumerate(TARGETS[modality])})
        return frame.sort_values(["target_order", "group_order", "donor_id"]).drop(
            columns=["group_order", "target_order"]
        )
    finally:
        adata.file.close()


def exact_label_permutation(pbc: np.ndarray, normal: np.ndarray) -> tuple[float, float]:
    combined = np.concatenate([pbc, normal])
    observed = float(pbc.mean() - normal.mean())
    effects = []
    for selected in combinations(range(len(combined)), len(pbc)):
        selected_mask = np.zeros(len(combined), dtype=bool)
        selected_mask[list(selected)] = True
        effects.append(float(combined[selected_mask].mean() - combined[~selected_mask].mean()))
    effects_array = np.asarray(effects)
    tolerance = 1e-12
    one_sided = float(np.mean(effects_array >= observed - tolerance))
    two_sided = float(np.mean(np.abs(effects_array) >= abs(observed) - tolerance))
    return one_sided, two_sided


def test_modality(frame: pd.DataFrame, modality: str) -> pd.DataFrame:
    rows: list[dict[str, float | int | str | bool]] = []
    for target in TARGETS[modality]:
        target_frame = frame[frame["target"] == target]
        maximum_threshold = 50 if modality == "scRNA" and target in {"NK_all", "T_NK_union"} else 20
        thresholds = [1, 10, 20] + ([50] if maximum_threshold == 50 else [])
        for minimum_cells in thresholds:
            eligible = target_frame[target_frame["n_cells"] >= minimum_cells]
            pbc = eligible.loc[eligible["group"] == "PBC", "logCPM"].to_numpy(float)
            normal = eligible.loc[eligible["group"] == "Normal", "logCPM"].to_numpy(float)
            if len(pbc) < 2 or len(normal) < 2:
                continue
            effect = float(pbc.mean() - normal.mean())
            exact_one_sided, exact_two_sided = exact_label_permutation(pbc, normal)
            welch = stats.ttest_ind(pbc, normal, equal_var=False, alternative="greater")
            rows.append(
                {
                    "modality": modality,
                    "target": target,
                    "min_cells": minimum_cells,
                    "n_pbc": len(pbc),
                    "n_normal": len(normal),
                    "pbc_logCPM": ";".join(f"{value:.4f}" for value in pbc),
                    "normal_logCPM": ";".join(f"{value:.4f}" for value in normal),
                    "effect_logCPM": effect,
                    "exact_one_sided": exact_one_sided,
                    "exact_two_sided": exact_two_sided,
                    "welch_one_sided": float(welch.pvalue),
                    "pbc_all_above_max_normal": bool(pbc.min() > normal.max()),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    download_if_requested(args.scrna_h5ad, "scRNA", args.download)
    download_if_requested(args.snrna_h5ad, "snRNA", args.download)
    OUT.mkdir(parents=True, exist_ok=True)

    donor_frames = []
    test_frames = []
    for modality, path in [("scRNA", args.scrna_h5ad), ("snRNA", args.snrna_h5ad)]:
        donor = aggregate_modality(path, modality)
        tests = test_modality(donor, modality)
        donor.to_csv(OUT / f"{modality}_IL12RB2_raw_pseudobulk_donor.tsv", sep="\t", index=False)
        tests.to_csv(OUT / f"{modality}_IL12RB2_raw_pseudobulk_tests.tsv", sep="\t", index=False)
        donor_frames.append(donor)
        test_frames.append(tests)

    pd.concat(donor_frames, ignore_index=True).to_csv(
        OUT / "andrews_IL12RB2_raw_pseudobulk_donor_all_modalities.tsv", sep="\t", index=False
    )
    pd.concat(test_frames, ignore_index=True).to_csv(
        OUT / "andrews_IL12RB2_raw_pseudobulk_tests_all_modalities.tsv", sep="\t", index=False
    )
    print(f"Wrote Andrews donor-level outputs to {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
