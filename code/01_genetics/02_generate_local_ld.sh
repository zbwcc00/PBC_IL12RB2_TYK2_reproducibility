#!/usr/bin/env bash
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
plink2="${PLINK2_BIN:-plink2}"
audit="$repo/work/genetics/fine_mapping/input/pbc_locus_allele_alignment_audit.tsv"
out_dir="$repo/work/genetics/fine_mapping/ld"
manifest="$out_dir/pbc_locus_local_ld_manifest.tsv"

mkdir -p "$out_dir"
printf 'locus_id\tvariant_n\tld_prefix\tld_matrix\tld_vars\tstatus\n' > "$manifest"

tail -n +2 "$audit" | tr -d '\r' | while IFS=$'\t' read -r \
  locus_id priority pbc_locus lead_rsid panel_n overlap_n coordinate_n harmonized_n \
  alignment_counts data_path snplist_path reference_prefix status; do
  [[ "$status" == "ready" ]] || continue
  pfile="$reference_prefix"
  snplist="$snplist_path"
  output_prefix="$out_dir/$locus_id"

  "$plink2" \
    --pfile "$pfile" \
    --extract "$snplist" \
    --r-unphased square bin4 ref-based \
    --threads 8 \
    --out "$output_prefix"

  matrix="$output_prefix.unphased.vcor1.bin"
  vars="$output_prefix.unphased.vcor1.bin.vars"
  [[ -s "$matrix" && -s "$vars" ]]
  observed_n="$(wc -l < "$vars")"
  [[ "$observed_n" -eq "$harmonized_n" ]]
  printf '%s\t%s\t%s\t%s\t%s\tready\n' \
    "$locus_id" "$observed_n" "$output_prefix" "$matrix" "$vars" >> "$manifest"
done

printf 'COMPLETE\t%s\n' "$(($(wc -l < "$manifest") - 1))"
