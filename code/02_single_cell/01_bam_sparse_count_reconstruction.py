from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pysam
from scipy import sparse
from scipy.io import mmwrite


TARGETS = {"IL12RB2", "IL12RB1", "JAK1", "JAK2", "TYK2", "STAT4", "IFNG", "IL12A", "IL12B"}
UMI_COUNT_FLAG = 8


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--donor", required=True)
    parser.add_argument("--condition", choices=["PBC", "control"], required=True)
    parser.add_argument("--expected-md5", required=True)
    parser.add_argument("--skip-md5", action="store_true")
    parser.add_argument("--progress-every", type=int, default=25_000_000)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    if not args.skip_md5:
        observed_md5 = md5sum(args.bam)
        if observed_md5.lower() != args.expected_md5.lower():
            raise RuntimeError(f"MD5 mismatch: expected {args.expected_md5}, observed {observed_md5}")
    else:
        observed_md5 = "skipped"

    counts: Counter[tuple[str, str]] = Counter()
    genes: dict[str, str] = {}
    audit = Counter()
    barcodes_seen = set()
    target_read_audit = Counter()
    started = time.monotonic()
    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        header = bam.header.to_dict()
        for read in bam.fetch(until_eof=True):
            audit["all_alignments"] += 1
            if args.progress_every > 0 and audit["all_alignments"] % args.progress_every == 0:
                print(json.dumps({
                    "sample": args.sample,
                    "alignments_scanned": audit["all_alignments"],
                    "counted_UMIs": audit["counted_UMIs"],
                    "elapsed_seconds": round(time.monotonic() - started, 1),
                }), flush=True)
            if read.is_secondary:
                audit["secondary_alignments"] += 1
                continue
            if read.is_supplementary:
                audit["supplementary_alignments"] += 1
                continue
            try:
                barcode = read.get_tag("CB")
                umi = read.get_tag("UB")
            except KeyError:
                audit["primary_missing_CB_or_UB"] += 1
                continue
            barcodes_seen.add(barcode)
            try:
                flags = int(read.get_tag("xf"))
            except KeyError:
                audit["primary_missing_xf"] += 1
                continue
            if flags & UMI_COUNT_FLAG == 0:
                audit["primary_not_UMI_COUNT"] += 1
                continue
            try:
                gene_ids = str(read.get_tag("GX")).split(";")
                gene_names = str(read.get_tag("GN")).split(";")
            except KeyError:
                audit["UMI_COUNT_missing_GX_or_GN"] += 1
                continue
            if len(gene_ids) != 1 or len(gene_names) != 1:
                audit["UMI_COUNT_multigene"] += 1
                continue
            gene_id, gene_name = gene_ids[0], gene_names[0]
            genes.setdefault(gene_id, gene_name)
            counts[(barcode, gene_id)] += 1
            audit["counted_UMIs"] += 1
            if gene_name.upper() in TARGETS:
                target_read_audit[(gene_name.upper(), barcode)] += 1

    barcodes = sorted({barcode for barcode, _ in counts})
    gene_ids = sorted(genes)
    barcode_index = {barcode: index for index, barcode in enumerate(barcodes)}
    gene_index = {gene_id: index for index, gene_id in enumerate(gene_ids)}
    rows, columns, values = [], [], []
    for (barcode, gene_id), value in counts.items():
        rows.append(gene_index[gene_id])
        columns.append(barcode_index[barcode])
        values.append(value)
    matrix = sparse.coo_matrix((values, (rows, columns)), shape=(len(gene_ids), len(barcodes)), dtype="int64")
    with gzip.open(args.output / "matrix.mtx.gz", "wb") as handle:
        mmwrite(handle, matrix, field="integer")
    with gzip.open(args.output / "barcodes.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        for barcode in barcodes:
            handle.write(f"{barcode}\n")
    with gzip.open(args.output / "features.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        for gene_id in gene_ids:
            handle.write(f"{gene_id}\t{genes[gene_id]}\tGene Expression\n")

    barcode_totals = pd.DataFrame({
        "barcode": barcodes,
        "total_UMI": list(map(int, sparse.csc_matrix(matrix).sum(axis=0).A1)),
        "detected_genes": list(map(int, sparse.csc_matrix(matrix).getnnz(axis=0))),
    })
    barcode_totals.to_csv(args.output / "barcode_QC.tsv.gz", sep="\t", index=False, compression="gzip")
    target_rows = []
    for (gene, barcode), value in sorted(target_read_audit.items()):
        target_rows.append({"gene": gene, "barcode": barcode, "UMI": value})
    pd.DataFrame(target_rows, columns=["gene", "barcode", "UMI"]).to_csv(
        args.output / "target_gene_barcode_counts.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame([{"metric": key, "value": value} for key, value in sorted(audit.items())]).to_csv(
        args.output / "extraction_audit.tsv", sep="\t", index=False
    )

    pg = header.get("PG", [])
    manifest = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "bam": str(args.bam),
        "bam_bytes": args.bam.stat().st_size,
        "bam_md5": observed_md5,
        "expected_md5": args.expected_md5,
        "sample": args.sample,
        "donor": args.donor,
        "condition": args.condition,
        "bam_decompression_threads": args.threads,
        "count_rule": "primary, nonsupplementary alignment with CB, UB, a unique GX/GN, and Cell Ranger xf UMI_COUNT bit 8; each such representative alignment contributes one UMI",
        "cellranger_source_reference": "10XGenomics/cellranger cr_bam ExtraFlags UMI_COUNT=8; repository snapshot stored under project 02_reference resources",
        "n_barcodes_with_counts": len(barcodes),
        "n_features_with_counts": len(gene_ids),
        "matrix_total_UMI": int(matrix.sum()),
        "BAM_programs": pg,
        "raw_unfiltered_caveat": "CB is the corrected barcode tag. Barcodes without CB cannot be reconstructed into a standard corrected-barcode matrix from this rule and are retained only in audit metrics.",
    }
    (args.output / "extraction_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
