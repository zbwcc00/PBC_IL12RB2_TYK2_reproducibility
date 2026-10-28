from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "molecular_qtl" / "mesa"
OUT = REPO / "work" / "molecular_qtl" / "mesa"
OLINK = RAW / "olink_cis_trans_pQTL_combined.tsv.gz"
SOMASCAN = RAW / "soma_cis_trans_pQTL_combined.tsv.gz"
GENE = "IL12RB2"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce the MESA same-sample Olink-SomaScan IL12RB2 cis-pQTL comparison."
    )
    parser.add_argument("--olink", type=Path, default=OLINK)
    parser.add_argument("--somascan", type=Path, default=SOMASCAN)
    return parser.parse_args()


def read_target(path: Path, platform: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Download the MESA Olink/SomaScan release from "
            "https://doi.org/10.5281/zenodo.17642645 and place it under "
            "data/raw/molecular_qtl/mesa/."
        )
    rows: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if GENE in row.get("phenotype_id", "").upper() and row.get("cis_trans") == "cis":
                rows.append(row)
    if not rows:
        raise RuntimeError(f"No cis {GENE} records found in {path.name}")
    frame = pd.DataFrame(rows)
    numeric = [
        "chromosome",
        "base_pair_location",
        "beta",
        "standard_error",
        "effect_allele_frequency",
        "p_value",
    ]
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["platform"] = platform
    frame["variant_key"] = (
        frame["chromosome"].astype("Int64").astype(str)
        + ":"
        + frame["base_pair_location"].astype("Int64").astype(str)
    )
    return frame.sort_values("p_value").drop_duplicates("variant_key", keep="first")


def align_platforms(olink: pd.DataFrame, somascan: pd.DataFrame) -> pd.DataFrame:
    olink = olink.rename(
        columns={
            "beta": "beta_o",
            "standard_error": "se_o",
            "effect_allele": "ea_o",
            "other_allele": "oa_o",
            "effect_allele_frequency": "eaf_o",
            "p_value": "p_o",
            "variant_id": "variant_id_o",
            "phenotype_id": "phenotype_o",
        }
    )
    somascan = somascan.rename(
        columns={
            "beta": "beta_s",
            "standard_error": "se_s",
            "effect_allele": "ea_s",
            "other_allele": "oa_s",
            "effect_allele_frequency": "eaf_s",
            "p_value": "p_s",
            "variant_id": "variant_id_s",
            "phenotype_id": "phenotype_s",
        }
    )
    merged = olink.merge(somascan, on="variant_key", how="inner", suffixes=("", "_sraw"))

    def alignment(row: pd.Series) -> int:
        if row.ea_o == row.ea_s and row.oa_o == row.oa_s:
            return 1
        if row.ea_o == row.oa_s and row.oa_o == row.ea_s:
            return -1
        return 0

    merged["allele_alignment"] = merged.apply(alignment, axis=1)
    merged = merged[merged["allele_alignment"] != 0].copy()
    merged["beta_s_aligned"] = merged["beta_s"] * merged["allele_alignment"]
    merged["direction_concordant"] = np.sign(merged["beta_o"]) == np.sign(merged["beta_s_aligned"])
    merged["z_o"] = merged["beta_o"] / merged["se_o"]
    merged["z_s"] = merged["beta_s_aligned"] / merged["se_s"]
    merged["z_product"] = merged["z_o"] * merged["z_s"]
    merged["min_p"] = merged[["p_o", "p_s"]].min(axis=1)
    return merged.sort_values("min_p")


def coloc_abf(frame: pd.DataFrame, p1: float, p2: float, p12: float, prior_variance: float = 0.04) -> dict[str, float | int | str]:
    valid = frame[
        np.isfinite(frame["beta_o"])
        & np.isfinite(frame["se_o"])
        & np.isfinite(frame["beta_s_aligned"])
        & np.isfinite(frame["se_s"])
    ].copy()
    beta_1 = valid["beta_o"].to_numpy(float)
    beta_2 = valid["beta_s_aligned"].to_numpy(float)
    se_1 = valid["se_o"].to_numpy(float)
    se_2 = valid["se_s"].to_numpy(float)
    shrink_1 = prior_variance / (prior_variance + se_1**2)
    shrink_2 = prior_variance / (prior_variance + se_2**2)
    log_bf_1 = 0.5 * (np.log1p(-shrink_1) + shrink_1 * (beta_1 / se_1) ** 2)
    log_bf_2 = 0.5 * (np.log1p(-shrink_2) + shrink_2 * (beta_2 / se_2) ** 2)
    bf_1 = np.exp(np.clip(log_bf_1, -700, 700))
    bf_2 = np.exp(np.clip(log_bf_2, -700, 700))
    shared_sum = np.sum(bf_1 * bf_2)
    masses = np.array(
        [
            1.0,
            p1 * np.sum(bf_1),
            p2 * np.sum(bf_2),
            p1 * p2 * (np.sum(bf_1) * np.sum(bf_2) - shared_sum),
            p12 * shared_sum,
        ]
    )
    posterior = masses / masses.sum()
    top_index = int(np.argmax(bf_1 * bf_2))
    return {
        "p1": p1,
        "p2": p2,
        "p12": p12,
        "W": prior_variance,
        "n_variants": int(len(valid)),
        "PP.H0": float(posterior[0]),
        "PP.H1": float(posterior[1]),
        "PP.H2": float(posterior[2]),
        "PP.H3": float(posterior[3]),
        "PP.H4": float(posterior[4]),
        "top_shared_variant": str(valid.iloc[top_index]["variant_key"]),
        "top_shared_product_logBF": float(np.log(bf_1[top_index] * bf_2[top_index])),
    }


def main() -> None:
    args = parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    olink = read_target(args.olink, "Olink Explore 3072")
    somascan = read_target(args.somascan, "SomaScan 7k")
    aligned = align_platforms(olink, somascan)
    aligned_path = OUT / "28_MESA_IL12RB2_shared_cis_variants_aligned.tsv"
    aligned.to_csv(aligned_path, sep="\t", index=False)

    prior_grid = [
        (1e-4, 1e-4, 1e-5),
        (1e-4, 1e-4, 1e-6),
        (1e-5, 1e-5, 1e-6),
        (1e-3, 1e-3, 1e-5),
    ]
    coloc = pd.DataFrame([coloc_abf(aligned, *priors) for priors in prior_grid])
    coloc.to_csv(OUT / "31_MESA_IL12RB2_coloc_abf_single_causal_sensitivity.tsv", sep="\t", index=False)

    summary = {
        "olink_unique_variants": int(len(olink)),
        "somascan_unique_variants": int(len(somascan)),
        "allele_aligned_shared_variants": int(len(aligned)),
        "direction_concordant_n": int(aligned["direction_concordant"].sum()),
        "direction_concordant_fraction": float(aligned["direction_concordant"].mean()),
        "pearson_z": float(aligned[["z_o", "z_s"]].corr().iloc[0, 1]),
        "spearman_z": float(aligned[["z_o", "z_s"]].corr(method="spearman").iloc[0, 1]),
        "maximum_single_causal_PP.H4": float(coloc["PP.H4"].max()),
        "interpretation": "Exploratory single-causal ABF sensitivity analysis; not coloc.susie.",
    }
    (OUT / "29_MESA_IL12RB2_cross_platform_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
