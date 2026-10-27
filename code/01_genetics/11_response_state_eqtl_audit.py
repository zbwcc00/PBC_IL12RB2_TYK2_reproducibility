from __future__ import annotations

import csv
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
import requests


REPO = Path(__file__).resolve().parents[2]
WORK_DIR = REPO / "work" / "genetics" / "response_state_eqtl"
RAW_DIR = WORK_DIR / "raw" / "permuted"
PROCESSED_DIR = WORK_DIR / "processed"
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

METADATA_PATH = REPO / "data" / "raw" / "external_eqtl" / "dataset_metadata_r8_beta.tsv"
GENE_ID = "ENSG00000081985"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path):
    if path.exists() and path.stat().st_size > 0:
        return "cached"
    temporary = path.with_suffix(path.suffix + ".part")
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for block in response.iter_content(chunk_size=1024 * 1024):
                if block:
                    handle.write(block)
    temporary.replace(path)
    return "downloaded"


def find_gene_column(columns: list[str]) -> str:
    for candidate in (
        "molecular_trait_id",
        "gene_id",
        "phenotype_id",
        "molecular_trait_object_id",
    ):
        if candidate in columns:
            return candidate
    raise RuntimeError(f"No gene identifier column in {columns}")


def first_present(row: pd.Series, candidates: tuple[str, ...]):
    for candidate in candidates:
        if candidate in row.index:
            return row[candidate]
    return ""


def write_tsv(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    metadata = pd.read_csv(METADATA_PATH, sep="\t", dtype=str).fillna("")
    randolph = metadata[
        (metadata["study_id"] == "QTS000036")
        & metadata["sample_group"].isin(
            [
                "B_NI",
                "B_flu",
                "CD4_T_NI",
                "CD4_T_flu",
                "CD8_T_NI",
                "CD8_T_flu",
                "NK_NI",
                "NK_flu",
            ]
        )
    ]
    cytoimmgen = metadata[metadata["study_id"] == "QTS000041"]
    selected = pd.concat([randolph, cytoimmgen], ignore_index=True)

    con = duckdb.connect()
    audit_rows = []
    schema_by_study = {}
    total = len(selected)
    for index, meta in selected.iterrows():
        study_id = meta["study_id"]
        dataset_id = meta["dataset_id"]
        url = (
            "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/r8_beta/sumstats/"
            f"{study_id}/{dataset_id}/{dataset_id}.permuted.parquet"
        )
        local_path = RAW_DIR / f"{dataset_id}.permuted.parquet"
        download_status = download(url, local_path)
        schema = con.execute(
            "DESCRIBE SELECT * FROM read_parquet(?)", [str(local_path)]
        ).fetchdf()
        columns = schema["column_name"].tolist()
        schema_by_study.setdefault(study_id, columns)
        gene_column = find_gene_column(columns)
        gene_rows = con.execute(
            f"SELECT * FROM read_parquet(?) WHERE starts_with({gene_column}, ?)",
            [str(local_path), GENE_ID],
        ).fetchdf()

        base = {
            "study_id": study_id,
            "dataset_id": dataset_id,
            "study_label": meta["study_label"],
            "sample_group": meta["sample_group"],
            "tissue_label": meta["tissue_label"],
            "condition_label": meta["condition_label"],
            "sample_size": meta["sample_size"],
            "pmid": meta["pmid"],
            "download_status": download_status,
            "file_bytes": local_path.stat().st_size,
            "sha256": sha256(local_path),
            "gene_rows": len(gene_rows),
        }
        if gene_rows.empty:
            audit_rows.append(
                {
                    **base,
                    "molecular_trait_id": "",
                    "top_variant": "",
                    "beta": "",
                    "se": "",
                    "pvalue": "",
                    "p_beta": "",
                    "p_perm": "",
                    "q_value": "",
                    "significant_p_beta_0_05": False,
                    "credible_signal_candidate": False,
                }
            )
        else:
            for _, row in gene_rows.iterrows():
                p_beta = first_present(row, ("p_beta", "pval_beta"))
                try:
                    p_beta_numeric = float(p_beta)
                except (TypeError, ValueError):
                    p_beta_numeric = float("nan")
                audit_rows.append(
                    {
                        **base,
                        "molecular_trait_id": row[gene_column],
                        "top_variant": first_present(
                            row,
                            (
                                "variant_id",
                                "variant",
                                "top_variant",
                                "snp_id",
                            ),
                        ),
                        "beta": first_present(row, ("beta", "slope")),
                        "se": first_present(row, ("se", "slope_se")),
                        "pvalue": first_present(
                            row, ("pvalue", "p_nominal", "pval_nominal")
                        ),
                        "p_beta": p_beta,
                        "p_perm": first_present(
                            row, ("p_perm", "pvalue_perm", "empirical_feature_p_value")
                        ),
                        "q_value": first_present(
                            row, ("q_value", "qvalue", "global_corrected_p_value")
                        ),
                        "significant_p_beta_0_05": bool(p_beta_numeric < 0.05),
                        "credible_signal_candidate": bool(p_beta_numeric < 0.05),
                    }
                )
        if (index + 1) % 10 == 0 or index + 1 == total:
            print(f"processed {index + 1}/{total}", flush=True)

    write_tsv(PROCESSED_DIR / "01_IL12RB2_response_state_gene_level_audit.tsv", audit_rows)
    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "gene_id": GENE_ID,
        "metadata_path": str(METADATA_PATH),
        "studies": {
            "QTS000036": {
                "label": "Randolph_2021",
                "contrast_context": "baseline versus influenza 6 h",
                "interpretation": "condition-specific eQTL, not automatically a formal genotype-by-stimulation interaction QTL",
            },
            "QTS000041": {
                "label": "Cytoimmgen",
                "contrast_context": "unstimulated and anti-CD3/CD28 CD4 T-cell states across time",
                "interpretation": "activation-state eQTL; IL-12 specificity cannot be inferred",
            },
        },
        "schemas": schema_by_study,
        "software": {"duckdb": duckdb.__version__, "pandas": pd.__version__},
        "elapsed_note": "Only gene-level permuted files were downloaded in this gate. Full regional files are queried only for nominated states.",
    }
    (PROCESSED_DIR / "02_gene_level_audit_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    started = time.time()
    main()
    print(f"elapsed_seconds={time.time() - started:.1f}", flush=True)
