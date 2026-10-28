from __future__ import annotations

import argparse
import math
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "single_cell" / "gse270739"
OUT = REPO / "work" / "single_cell" / "gse270739"
META_FILE = RAW / "figure2_input_metadata.csv"
UMI_FILE = RAW / "figure2_input_UMI.csv"
CLINICAL_FILE = RAW / "clinical_supplement.xlsx"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce the GSE270739 within-PBC antigen-context IL12RB2 analysis."
    )
    parser.add_argument("--metadata", type=Path, default=META_FILE)
    parser.add_argument("--expression", type=Path, default=UMI_FILE)
    parser.add_argument("--clinical", type=Path, default=CLINICAL_FILE)
    return parser.parse_args()


def patient_core(value: object) -> str:
    text = str(value).strip().replace("_", "-").replace(" ", "-")
    fields = [part for part in text.split("-") if part]
    return "-".join(fields[:2]) if len(fields) >= 2 else text


def exact_sign_flip_p(differences: pd.Series) -> float:
    values = differences.to_numpy(float)
    observed = abs(values.mean())
    permuted = [abs(np.mean(values * np.asarray(signs))) for signs in product([-1.0, 1.0], repeat=len(values))]
    return float(np.mean(np.asarray(permuted) >= observed - 1e-12))


def paired_summary(frame: pd.DataFrame, column: str) -> tuple[dict[str, float | int | str], pd.DataFrame, pd.Series]:
    wide = frame.pivot(index="patient_core", columns="condition", values=column).dropna()
    differences = wide["PDC-E2"] - wide["foreign antigen"]
    sample_size = len(differences)
    mean = float(differences.mean())
    standard_deviation = float(differences.std(ddof=1)) if sample_size > 1 else float("nan")
    standard_error = standard_deviation / math.sqrt(sample_size) if sample_size > 1 else float("nan")
    critical_value = stats.t.ppf(0.975, sample_size - 1) if sample_size > 1 else float("nan")
    paired_test = stats.ttest_rel(wide["PDC-E2"], wide["foreign antigen"])
    return {
        "metric": column,
        "n_independent_patients": sample_size,
        "mean_paired_difference": mean,
        "ci95_low": mean - critical_value * standard_error,
        "ci95_high": mean + critical_value * standard_error,
        "paired_t_p_two_sided": float(paired_test.pvalue),
        "exact_signflip_p_two_sided": exact_sign_flip_p(differences),
        "positive_direction_n": int((differences > 0).sum()),
        "zero_direction_n": int((differences == 0).sum()),
    }, wide, differences


def aggregate(frame: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, float | int | str]] = []
    for (patient, condition), group in frame.groupby(["patient_core", "condition"]):
        gene_total = float(group["IL12RB2"].sum())
        library_total = float(group["nCount_RNA"].sum())
        cpm = gene_total / library_total * 1e6 if library_total > 0 else np.nan
        records.append(
            {
                "patient_core": patient,
                "disease": group["disease"].iloc[0],
                "condition": condition,
                "n_cells": len(group),
                "IL12RB2_total": gene_total,
                "library_total": library_total,
                "IL12RB2_CPM": cpm,
                "IL12RB2_log2CPM": math.log2(cpm + 0.5),
                "IL12RB2_positive_fraction": float((group["IL12RB2"] > 0).mean()),
            }
        )
    return pd.DataFrame(records)


def require_inputs(paths: list[Path]) -> None:
    missing = [path for path in paths if not path.exists()]
    if missing:
        formatted = "\n".join(f"- {path}" for path in missing)
        raise FileNotFoundError(
            "Missing GSE270739 inputs:\n"
            f"{formatted}\n"
            "Place the processed expression files and clinical supplement under "
            "data/raw/single_cell/gse270739/. See data/metadata/data_source_registry.tsv."
        )


def main() -> None:
    args = parse_args()
    require_inputs([args.metadata, args.expression, args.clinical])
    OUT.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(args.metadata)
    expression = pd.read_csv(args.expression, usecols=["IL12RB2"])["IL12RB2"]
    if len(metadata) != len(expression):
        raise RuntimeError("Metadata and expression rows do not align")
    metadata = metadata.copy()
    metadata["IL12RB2"] = expression.to_numpy()
    metadata["patient_core"] = metadata["Patient_ID"].map(patient_core)

    clinical = pd.read_excel(args.clinical, header=2)
    clinical = clinical[clinical["sample"].notna()].copy()
    clinical["patient_core"] = clinical["sample"].map(patient_core)
    disease = clinical.drop_duplicates("patient_core").set_index("patient_core")["AILD"].to_dict()
    metadata["disease"] = metadata["patient_core"].map(disease)
    metadata["condition"] = np.where(
        metadata["Ag_reactivity"].astype(str).str.upper().str.startswith("PDCE2"),
        "PDC-E2",
        "foreign antigen",
    )

    audit = (
        metadata.groupby(["patient_core", "Patient_ID", "disease", "Ag_reactivity"], dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    audit.to_csv(OUT / "01_patient_antigen_cell_audit.tsv", sep="\t", index=False)

    primary = aggregate(metadata[metadata["disease"].astype(str).str.strip().eq("PBC")])
    eligible = primary.groupby("patient_core")["condition"].nunique()
    primary = primary[primary["patient_core"].isin(eligible[eligible == 2].index)].copy()
    primary.to_csv(OUT / "02_pure_PBC_patient_condition_pseudobulk.tsv", sep="\t", index=False)

    sensitivity = aggregate(
        metadata[metadata["disease"].astype(str).str.strip().isin(["PBC", "AIH/PBC"])]
    )
    eligible_sensitivity = sensitivity.groupby("patient_core")["condition"].nunique()
    sensitivity = sensitivity[
        sensitivity["patient_core"].isin(eligible_sensitivity[eligible_sensitivity == 2].index)
    ].copy()
    sensitivity.to_csv(OUT / "03_PBC_plus_overlap_sensitivity.tsv", sep="\t", index=False)

    summaries: list[dict[str, float | int | str]] = []
    primary_differences: pd.Series | None = None
    for analysis, frame in [("pure_PBC", primary), ("PBC_plus_AIH_PBC_overlap", sensitivity)]:
        for metric in ["IL12RB2_log2CPM", "IL12RB2_positive_fraction"]:
            summary, _, differences = paired_summary(frame, metric)
            summary["analysis"] = analysis
            summaries.append(summary)
            if analysis == "pure_PBC" and metric == "IL12RB2_log2CPM":
                primary_differences = differences
    pd.DataFrame(summaries).to_csv(OUT / "04_paired_inference_summary.tsv", sep="\t", index=False)

    if primary_differences is None:
        raise RuntimeError("Primary paired IL12RB2 analysis was not generated")
    leave_one_out = []
    for omitted in primary_differences.index:
        retained = primary_differences.drop(index=omitted)
        leave_one_out.append(
            {
                "omitted_patient": omitted,
                "n_retained": len(retained),
                "mean_paired_difference": float(retained.mean()),
                "direction_positive": bool(retained.mean() > 0),
            }
        )
    pd.DataFrame(leave_one_out).to_csv(OUT / "05_leave_one_patient_out.tsv", sep="\t", index=False)
    print(f"Wrote GSE270739 donor-level outputs to {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
