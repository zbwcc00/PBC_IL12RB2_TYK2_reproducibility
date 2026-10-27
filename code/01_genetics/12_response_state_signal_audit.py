from __future__ import annotations

import csv
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
import requests


REPO = Path(__file__).resolve().parents[2]
WORK_DIR = REPO / "work" / "genetics" / "response_state_eqtl"
RAW_DIR = WORK_DIR / "raw" / "susie"
PROCESSED_DIR = WORK_DIR / "processed"
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

METADATA_PATH = REPO / "data" / "raw" / "external_eqtl" / "dataset_metadata_r8_beta.tsv"
GENE_AUDIT_PATH = PROCESSED_DIR / "01_IL12RB2_response_state_gene_level_audit.tsv"
PBC_CS_PATH = (
    REPO
    / "work"
    / "genetics"
    / "regulatory_annotation"
    / "processed"
    / "01_credible_set_variant_regulatory_annotation.tsv"
)
GENE_ID = "ENSG00000081985"
STUDY_ID = "QTS000041"

PAIRINGS = {
    "QTD000691": ("QTD000694", "memory lineage; 40 h stimulated versus 16 h unstimulated"),
    "QTD000692": ("QTD000694", "memory lineage; 5 d stimulated versus 16 h unstimulated"),
    "QTD000701": ("QTD000734", "TEM lineage; 5 d stimulated versus 16 h unstimulated"),
    "QTD000708": ("QTD000707", "TCM lineage; 40 h stimulated versus 16 h unstimulated"),
    "QTD000712": ("QTD000689", "naive lineage proxy; TN_A 16 h stimulated versus aggregate naive 16 h unstimulated"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> str:
    if path.exists() and path.stat().st_size > 0:
        return "cached"
    temporary = path.with_suffix(path.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for block in response.iter_content(chunk_size=1024 * 1024):
                if block:
                    handle.write(block)
    temporary.replace(path)
    return "downloaded"


def write_tsv(path: Path, frame: pd.DataFrame):
    frame.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def query_all_parquet(
    dataset_id: str,
    positions: list[int],
    metadata_row: dict,
) -> tuple[str, pd.DataFrame, str]:
    url = (
        "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/r8_beta/sumstats/"
        f"{STUDY_ID}/{dataset_id}/{dataset_id}.all.parquet"
    )
    placeholders = ",".join("?" for _ in positions)
    query = f"""
        SELECT molecular_trait_id, chromosome, position, ref, alt, variant,
               rsid, maf, pvalue, beta, se, type
        FROM read_parquet(?)
        WHERE starts_with(molecular_trait_id, ?)
          AND chromosome = ?
          AND position IN ({placeholders})
        ORDER BY position, pvalue
    """
    last_error = ""
    for attempt in range(1, 4):
        connection = duckdb.connect()
        connection.execute("SET enable_progress_bar=false")
        try:
            frame = connection.execute(
                query,
                [url, GENE_ID, "1", *positions],
            ).fetchdf()
            for key, value in metadata_row.items():
                frame[key] = value
            frame["source_url"] = url
            frame["query_attempt"] = attempt
            return dataset_id, frame, ""
        except Exception as error:
            last_error = repr(error)
            time.sleep(attempt * 3)
        finally:
            connection.close()
    return dataset_id, pd.DataFrame(), last_error


def main():
    metadata = pd.read_csv(METADATA_PATH, sep="\t", dtype=str).fillna("")
    metadata = metadata[metadata["study_id"] == STUDY_ID].copy()
    metadata_by_id = metadata.set_index("dataset_id").to_dict("index")
    gene_audit = pd.read_csv(GENE_AUDIT_PATH, sep="\t", dtype=str).fillna("")
    gene_audit["p_beta_numeric"] = pd.to_numeric(gene_audit["p_beta"], errors="coerce")
    nominated = sorted(
        gene_audit.loc[gene_audit["p_beta_numeric"] < 0.05, "dataset_id"].tolist()
    )
    expected = sorted(PAIRINGS)
    if nominated != expected:
        raise RuntimeError(f"Unexpected nominated states: {nominated}; expected {expected}")

    pbc_cs = pd.read_csv(PBC_CS_PATH, sep="\t", dtype=str).fillna("")
    pbc_positions = sorted(pbc_cs["ncbi_grch38_pos1"].astype(int).unique().tolist())
    pbc_lookup = pbc_cs.set_index("ncbi_grch38_pos1").to_dict("index")

    all_dataset_ids = sorted(set(nominated) | {pair[0] for pair in PAIRINGS.values()})
    susie_rows = []
    susie_summary = []
    download_manifest = []
    for dataset_id in all_dataset_ids:
        url = (
            "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/r8_beta/susie/"
            f"{STUDY_ID}/{dataset_id}/{dataset_id}.credible_sets.parquet"
        )
        local_path = RAW_DIR / f"{dataset_id}.credible_sets.parquet"
        status = download(url, local_path)
        download_manifest.append(
            {
                "dataset_id": dataset_id,
                "url": url,
                "local_path": str(local_path),
                "download_status": status,
                "file_bytes": local_path.stat().st_size,
                "sha256": sha256(local_path),
            }
        )
        connection = duckdb.connect()
        frame = connection.execute(
            "SELECT * FROM read_parquet(?) WHERE starts_with(molecular_trait_id, ?)",
            [str(local_path), GENE_ID],
        ).fetchdf()
        connection.close()
        meta = metadata_by_id[dataset_id]
        if not frame.empty:
            frame.insert(0, "dataset_id", dataset_id)
            frame.insert(1, "sample_group", meta["sample_group"])
            frame.insert(2, "condition_label", meta["condition_label"])
            frame.insert(3, "sample_size", meta["sample_size"])
            frame["pbc_credible_variant"] = frame["position"].astype(str).isin(pbc_lookup)
            frame["pbc_component"] = frame["position"].astype(str).map(
                lambda value: pbc_lookup.get(value, {}).get("component", "")
            )
            frame["pbc_rsid"] = frame["position"].astype(str).map(
                lambda value: pbc_lookup.get(value, {}).get("rsid", "")
            )
            susie_rows.append(frame)
        overlap = frame[frame["position"].isin(pbc_positions)] if not frame.empty else frame
        susie_summary.append(
            {
                "dataset_id": dataset_id,
                "sample_group": meta["sample_group"],
                "condition_label": meta["condition_label"],
                "sample_size": meta["sample_size"],
                "gene_credible_set_rows": len(frame),
                "gene_credible_set_count": frame["cs_id"].nunique() if not frame.empty else 0,
                "pbc_credible_variant_overlap_count": len(overlap),
                "pbc_overlap_rsids": ";".join(sorted(overlap["pbc_rsid"].dropna().astype(str).unique()))
                if not overlap.empty
                else "",
                "max_pip_overall": frame["pip"].max() if not frame.empty else "",
                "max_pip_pbc_overlap": overlap["pip"].max() if not overlap.empty else "",
            }
        )
        print(f"susie {dataset_id}: rows={len(frame)}, overlap={len(overlap)}", flush=True)

    write_tsv(PROCESSED_DIR / "03_susie_download_manifest.tsv", pd.DataFrame(download_manifest))
    combined_susie = pd.concat(susie_rows, ignore_index=True) if susie_rows else pd.DataFrame()
    write_tsv(PROCESSED_DIR / "04_IL12RB2_response_state_susie_rows.tsv", combined_susie)
    write_tsv(PROCESSED_DIR / "05_IL12RB2_response_state_susie_summary.tsv", pd.DataFrame(susie_summary))

    remote_frames = []
    remote_errors = []
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {}
        for dataset_id in all_dataset_ids:
            meta = metadata_by_id[dataset_id]
            compact_meta = {
                "dataset_id": dataset_id,
                "sample_group": meta["sample_group"],
                "condition_label": meta["condition_label"],
                "sample_size": meta["sample_size"],
            }
            future = executor.submit(
                query_all_parquet,
                dataset_id,
                pbc_positions,
                compact_meta,
            )
            futures[future] = dataset_id
        for future in as_completed(futures):
            dataset_id, frame, error = future.result()
            if error:
                remote_errors.append({"dataset_id": dataset_id, "error": error})
                print(f"all {dataset_id}: ERROR", flush=True)
            else:
                remote_frames.append(frame)
                print(f"all {dataset_id}: rows={len(frame)}", flush=True)

    all_variants = pd.concat(remote_frames, ignore_index=True) if remote_frames else pd.DataFrame()
    if not all_variants.empty:
        all_variants["pbc_rsid"] = all_variants["position"].astype(str).map(
            lambda value: pbc_lookup.get(value, {}).get("rsid", "")
        )
        all_variants["pbc_component"] = all_variants["position"].astype(str).map(
            lambda value: pbc_lookup.get(value, {}).get("component", "")
        )
        all_variants["pbc_pip"] = all_variants["position"].astype(str).map(
            lambda value: pbc_lookup.get(value, {}).get("pbc_pip", "")
        )
    write_tsv(PROCESSED_DIR / "06_IL12RB2_pbc_variant_state_effects.tsv", all_variants)
    write_tsv(PROCESSED_DIR / "07_remote_query_errors.tsv", pd.DataFrame(remote_errors))

    comparison_rows = []
    for stimulated_id, (baseline_id, pairing_note) in PAIRINGS.items():
        stimulated = all_variants[all_variants["dataset_id"] == stimulated_id]
        baseline = all_variants[all_variants["dataset_id"] == baseline_id]
        merged = stimulated.merge(
            baseline,
            on=["position", "ref", "alt", "pbc_rsid", "pbc_component", "pbc_pip"],
            how="outer",
            suffixes=("_stimulated", "_baseline"),
        )
        for _, row in merged.iterrows():
            stimulated_beta = pd.to_numeric(row.get("beta_stimulated"), errors="coerce")
            baseline_beta = pd.to_numeric(row.get("beta_baseline"), errors="coerce")
            comparison_rows.append(
                {
                    "stimulated_dataset_id": stimulated_id,
                    "baseline_dataset_id": baseline_id,
                    "stimulated_sample_group": metadata_by_id[stimulated_id]["sample_group"],
                    "baseline_sample_group": metadata_by_id[baseline_id]["sample_group"],
                    "pairing_note": pairing_note,
                    "pbc_rsid": row.get("pbc_rsid", ""),
                    "pbc_component": row.get("pbc_component", ""),
                    "pbc_pip": row.get("pbc_pip", ""),
                    "position": row.get("position", ""),
                    "ref": row.get("ref", ""),
                    "alt": row.get("alt", ""),
                    "beta_stimulated": row.get("beta_stimulated", ""),
                    "se_stimulated": row.get("se_stimulated", ""),
                    "pvalue_stimulated": row.get("pvalue_stimulated", ""),
                    "beta_baseline": row.get("beta_baseline", ""),
                    "se_baseline": row.get("se_baseline", ""),
                    "pvalue_baseline": row.get("pvalue_baseline", ""),
                    "delta_beta_descriptive": stimulated_beta - baseline_beta
                    if pd.notna(stimulated_beta) and pd.notna(baseline_beta)
                    else "",
                    "same_sign": bool(stimulated_beta * baseline_beta > 0)
                    if pd.notna(stimulated_beta) and pd.notna(baseline_beta)
                    else "",
                }
            )
    write_tsv(
        PROCESSED_DIR / "08_stimulated_unstimulated_variant_comparison.tsv",
        pd.DataFrame(comparison_rows),
    )

    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "gene_id": GENE_ID,
        "study_id": STUDY_ID,
        "nominated_states": nominated,
        "comparison_pairings": PAIRINGS,
        "pbc_credible_positions_grch38": pbc_positions,
        "interpretation_limits": [
            "The nominated datasets are condition-specific eQTL maps, not formal genotype-by-stimulation interaction tests.",
            "Several stimulated-versus-unstimulated pairings differ in collection time or cell-state granularity.",
            "Beta differences are descriptive because covariance-aware interaction statistics are unavailable.",
            "The stimulation is anti-CD3/CD28 rather than IL-12-specific.",
        ],
        "software": {"duckdb": duckdb.__version__, "pandas": pd.__version__},
    }
    (PROCESSED_DIR / "09_signal_audit_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    started = time.time()
    main()
    print(f"elapsed_seconds={time.time() - started:.1f}", flush=True)
