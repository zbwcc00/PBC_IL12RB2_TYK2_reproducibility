from __future__ import annotations

import csv
import hashlib
import os
import re
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import duckdb


REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "external_eqtl"
WORK = REPO / "work" / "genetics" / "external_eqtl"
AUDIT = RAW / "audit"
CACHE = WORK / "cache"
OUT = WORK / "logbf"
LISTS = WORK / "manifests"
RESULTS = WORK / "results"
DUCKDB_HOME = CACHE / "duckdb_home"

TARGETS = [
    ("QTD000597", "FCRL3"),
    ("QTD000597", "FCRL3-AS1"),
    ("QTD000834", "FCRL3"),
    ("QTD000834", "FCRL3-AS1"),
    ("QTD000356", "FCRL3"),
    ("QTD000664", "FCRL3"),
    ("QTD000598", "IL12RB2"),
    ("QTD000356", "IL12RB2"),
    ("QTD000683", "IL12RB2"),
    ("QTD000356", "CD226"),
    ("QTD000683", "CD226"),
    ("QTD000781", "CD226"),
    ("QTD000689", "CD226"),
    ("QTD000326", "DOK6"),
    ("QTD000600", "DOK6"),
    ("QTD000683", "DOK6"),
    ("QTD000781", "DOK6"),
    ("QTD000694", "DOK6"),
    ("QTD000587", "DOK6"),
]

GENE_TO_LOCUS = {
    "FCRL3": "1q23.1",
    "FCRL3-AS1": "1q23.1",
    "IL12RB2": "1p31.3",
    "CD226": "18q22.2",
    "DOK6": "18q22.2",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sql_quote(value: str) -> str:
    return value.replace("'", "''")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repository_relative(path: Path) -> str:
    return path.resolve().relative_to(REPO).as_posix()


def component_number(cs_id: str) -> int:
    match = re.search(r"_L(\d+)$", cs_id)
    if not match:
        raise ValueError(f"Cannot parse SuSiE component from {cs_id}")
    return int(match.group(1))


def new_connection() -> duckdb.DuckDBPyConnection:
    connection = duckdb.connect(config={"home_directory": str(DUCKDB_HOME)})
    connection.execute("SET threads=1")
    connection.execute("SET http_timeout=300000")
    connection.execute("SET http_retries=5")
    return connection


def credible_rows(path: Path, gene_id: str) -> list[dict]:
    connection = new_connection()
    try:
        query = """
            SELECT molecular_trait_id, gene_id, cs_id, variant, rsid, pip, beta, se, z, region
            FROM read_parquet(?)
            WHERE molecular_trait_id = ?
            ORDER BY cs_id, pip DESC, variant
        """
        cursor = connection.execute(query, [str(path), gene_id])
        columns = [item[0] for item in cursor.description]
        return [dict(zip(columns, row)) for row in cursor.fetchall()]
    finally:
        connection.close()


def extract_dataset(dataset_id: str, genes: list[dict]) -> dict:
    url = genes[0]["lbf_variable_url"]
    gene_ids = sorted({row["ensembl_gene_id"] for row in genes})
    output = OUT / f"{dataset_id}__selected_lbf.tsv"
    conditions = []
    for row in genes:
        chromosome, interval = row["external_region"].removeprefix("chr").split(":", 1)
        start, end = map(int, interval.split("-", 1))
        conditions.append(
            "(molecular_trait_id = '{gene_id}' AND chromosome = '{chromosome}' "
            "AND position BETWEEN {start} AND {end})".format(
                gene_id=sql_quote(row["ensembl_gene_id"]),
                chromosome=sql_quote(chromosome),
                start=start,
                end=end,
            )
        )
    component_columns = sorted({
        f"lbf_variable{component}"
        for row in genes
        for component in row["external_cs_components"].split(";")
    })
    columns = [
        "molecular_trait_id",
        "region",
        "variant",
        *component_columns,
    ]
    query = (
        f"SELECT {','.join(columns)} FROM read_parquet('{sql_quote(url)}') "
        f"WHERE {' OR '.join(conditions)} ORDER BY molecular_trait_id, position, variant"
    )
    copy_sql = (
        f"COPY ({query}) TO '{sql_quote(output.as_posix())}' "
        "(FORMAT CSV, HEADER TRUE, DELIMITER '\t')"
    )
    last_error = None
    for attempt in range(1, 5):
        try:
            if output.exists():
                output.unlink()
            connection = new_connection()
            try:
                connection.execute(copy_sql)
                counts = dict(
                    connection.execute(
                        "SELECT molecular_trait_id, count(*) FROM read_csv_auto(?, delim='\\t', header=true) GROUP BY 1",
                        [str(output)],
                    ).fetchall()
                )
            finally:
                connection.close()
            missing = [gene_id for gene_id in gene_ids if counts.get(gene_id, 0) == 0]
            if missing:
                raise RuntimeError(f"No log-BF rows for {dataset_id}: {missing}")
            return {
                "dataset_id": dataset_id,
                "status": "success",
                "attempt": attempt,
                "path": repository_relative(output),
                "sha256": sha256(output),
                "bytes": output.stat().st_size,
                "source_path": url,
                "source_bytes": 0,
                "source_sha256": "remote_range_read_not_full_download",
                "counts": counts,
            }
        except Exception as error:
            last_error = repr(error)
            time.sleep(min(10 * attempt, 30))
    return {
        "dataset_id": dataset_id,
        "status": "failed",
        "attempt": 4,
        "path": repository_relative(output),
        "sha256": "",
        "bytes": 0,
        "source_path": url,
        "source_bytes": 0,
        "source_sha256": "remote_range_read_failed",
        "counts": {},
        "error": last_error,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    LISTS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    DUCKDB_HOME.mkdir(parents=True, exist_ok=True)

    coverage = read_tsv(AUDIT / "candidate_gene_external_coverage.tsv")
    downloads = read_tsv(AUDIT / "downloadable_dataset_manifest.tsv")
    coverage_index = {(row["dataset_id"], row["gene_symbol"]): row for row in coverage}
    download_index = {row["resource_id"]: row for row in downloads if row["database"] == "eQTL Catalogue"}

    dataset_filter = {
        value for value in os.getenv("DATASET_FILTER", "").split(",") if value
    }
    target_pairs = [
        pair for pair in TARGETS if not dataset_filter or pair[0] in dataset_filter
    ]
    selected = []
    cs_detail = []
    for order, key in enumerate(target_pairs, start=1):
        if key not in coverage_index:
            raise KeyError(f"Missing audited target: {key}")
        coverage_row = coverage_index[key]
        if int(coverage_row["susie_credible_set_count"] or 0) < 1:
            raise ValueError(f"Target lacks external credible set: {key}")
        dataset_id, gene_symbol = key
        download_row = download_index[dataset_id]
        if not download_row["lbf_variable_url"].endswith(".parquet"):
            raise ValueError(f"Primary signal-level extraction requires release 8 parquet: {dataset_id}")
        credible_path = CACHE / "susie_files" / f"{dataset_id}.credible_sets.parquet"
        if not credible_path.exists():
            raise FileNotFoundError(credible_path)
        gene_id = coverage_row["ensembl_gene_id"]
        rows = credible_rows(credible_path, gene_id)
        if not rows:
            raise ValueError(f"No credible-set rows for {dataset_id} {gene_symbol}")
        components = sorted({component_number(str(row["cs_id"])) for row in rows})
        regions = sorted({str(row["region"]) for row in rows})
        if len(regions) != 1:
            raise ValueError(f"Expected one cis region for {dataset_id} {gene_symbol}: {regions}")
        selected_row = {
            "analysis_order": order,
            "locus": GENE_TO_LOCUS[gene_symbol],
            "dataset_id": dataset_id,
            "study_id": coverage_row["study_id"],
            "study_label": coverage_row["study_label"],
            "sample_group": coverage_row["sample_group"],
            "sample_size": coverage_row["sample_size"],
            "ancestry": coverage_row["ancestry"],
            "release": coverage_row["release"],
            "genome_build": coverage_row["genome_build"],
            "gene_symbol": gene_symbol,
            "ensembl_gene_id": gene_id,
            "external_cs_components": ";".join(map(str, components)),
            "external_cs_count": len(components),
            "external_region": regions[0],
            "credible_sets_path": repository_relative(credible_path),
            "credible_sets_sha256": sha256(credible_path),
            "lbf_variable_url": download_row["lbf_variable_url"],
            "license": download_row["license"],
            "selection_rule": "one representative tissue/cell per independent study where possible",
        }
        selected.append(selected_row)
        for row in rows:
            cs_detail.append(
                {
                    **{name: selected_row[name] for name in [
                        "locus", "dataset_id", "study_id", "study_label", "sample_group",
                        "sample_size", "ancestry", "gene_symbol", "ensembl_gene_id"
                    ]},
                    "component": component_number(str(row["cs_id"])),
                    **row,
                }
            )

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in selected:
        grouped[row["dataset_id"]].append(row)

    dataset_results = {}
    max_workers = int(os.getenv("MAX_WORKERS", "3"))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(extract_dataset, dataset_id, rows): dataset_id
            for dataset_id, rows in grouped.items()
        }
        for future in as_completed(futures):
            result = future.result()
            dataset_results[result["dataset_id"]] = result
            print(
                f"{result['dataset_id']}\t{result['status']}\t"
                f"attempt={result['attempt']}\tbytes={result['bytes']}",
                flush=True,
            )

    failures = [result for result in dataset_results.values() if result["status"] != "success"]
    for row in selected:
        result = dataset_results[row["dataset_id"]]
        row["lbf_selected_path"] = result["path"]
        row["lbf_selected_sha256"] = result["sha256"]
        row["lbf_selected_bytes"] = result["bytes"]
        row["lbf_source_path"] = result["source_path"]
        row["lbf_source_bytes"] = result["source_bytes"]
        row["lbf_source_sha256"] = result["source_sha256"]
        row["lbf_gene_row_count"] = result["counts"].get(row["ensembl_gene_id"], 0)
        row["extraction_status"] = result["status"]
        row["extraction_error"] = result.get("error", "")

    write_tsv(
        LISTS / "external_signal_level_manifest.tsv",
        selected,
        list(selected[0].keys()),
    )
    write_tsv(
        RESULTS / "external_eqtl_susie_credible_set_details.tsv",
        cs_detail,
        list(cs_detail[0].keys()),
    )
    write_tsv(
        LISTS / "external_logbf_extraction_audit.tsv",
        [
            {
                "dataset_id": key,
                "status": value["status"],
                "attempt": value["attempt"],
                "path": value["path"],
                "sha256": value["sha256"],
                "bytes": value["bytes"],
                "source_path": value["source_path"],
                "source_bytes": value["source_bytes"],
                "source_sha256": value["source_sha256"],
                "error": value.get("error", ""),
            }
            for key, value in sorted(dataset_results.items())
        ],
        [
            "dataset_id", "status", "attempt", "path", "sha256", "bytes",
            "source_path", "source_bytes", "source_sha256", "error"
        ],
    )
    if failures:
        raise RuntimeError(f"External log-BF extraction failed for {len(failures)} datasets")


if __name__ == "__main__":
    main()
