from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "genetics"
WORK = REPO / "work" / "genetics" / "fine_mapping"
PANEL_MANIFEST = RAW / "1000_genomes" / "eur_signal_aware_panel_manifest.tsv"
PBC_GWAS = RAW / "GCST90061440_buildGRCh37.tsv"
OUT_DIR = WORK / "input"
AUDIT_OUT = OUT_DIR / "pbc_locus_allele_alignment_audit.tsv"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def resolve_manifest_path(value: str) -> Path:
    path = Path(str(value).strip())
    return path if path.is_absolute() else REPO / path


def align_to_alt(effect: str, other: str, reference: str, alternate: str) -> tuple[int | None, str]:
    effect = str(effect).upper()
    other = str(other).upper()
    reference = str(reference).upper()
    alternate = str(alternate).upper()
    if effect == alternate and other == reference:
        return 1, "direct_effect_alt"
    if effect == reference and other == alternate:
        return -1, "direct_effect_ref"
    if {effect, other} in ({"A", "T"}, {"C", "G"}):
        return None, "palindromic_non_direct_excluded"
    effect_comp = effect.translate(COMPLEMENT)
    other_comp = other.translate(COMPLEMENT)
    if effect_comp == alternate and other_comp == reference:
        return 1, "strand_complement_effect_alt"
    if effect_comp == reference and other_comp == alternate:
        return -1, "strand_complement_effect_ref"
    return None, "allele_mismatch"


def load_reference(prefix: Path) -> pd.DataFrame:
    pvar = pd.read_csv(Path(f"{prefix}.pvar"), sep="\t", dtype=str)
    pvar = pvar.rename(
        columns={
            "#CHROM": "chromosome_grch37",
            "POS": "position_grch37_1based",
            "ID": "rsid",
            "REF": "reference_allele",
            "ALT": "alternate_allele",
        }
    )
    pvar["position_grch37_1based"] = pd.to_numeric(pvar["position_grch37_1based"])
    afreq = pd.read_csv(Path(f"{prefix}.afreq"), sep="\t", dtype=str)
    afreq = afreq.rename(columns={"ID": "rsid", "ALT_FREQS": "reference_alt_freq"})
    afreq["reference_alt_freq"] = pd.to_numeric(afreq["reference_alt_freq"])
    return pvar.merge(afreq[["rsid", "reference_alt_freq"]], on="rsid", how="left")


def load_pbc(requested_rsids: set[str]) -> pd.DataFrame:
    usecols = [
        "variant_id",
        "chromosome",
        "base_pair_location",
        "effect_allele",
        "other_allele",
        "p_value",
        "beta",
        "standard_error",
    ]
    frames = []
    for chunk in pd.read_csv(PBC_GWAS, sep="\t", usecols=usecols, chunksize=750_000):
        matched = chunk.loc[chunk["variant_id"].isin(requested_rsids)]
        if not matched.empty:
            frames.append(matched.copy())
    pbc = pd.concat(frames, ignore_index=True)
    pbc = pbc.loc[pbc["beta"].notna() & pbc["standard_error"].gt(0)].copy()
    return pbc.sort_values(["p_value", "variant_id"]).drop_duplicates("variant_id", keep="first")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(PANEL_MANIFEST, sep="\t")
    references: dict[str, pd.DataFrame] = {}
    prefixes: dict[str, Path] = {}
    requested_rsids: set[str] = set()
    for row in manifest.to_dict("records"):
        locus = str(row["pbc_locus"])
        prefix = resolve_manifest_path(row["pgen_prefix"])
        reference = load_reference(prefix)
        prefixes[locus] = prefix
        references[locus] = reference
        requested_rsids.update(reference["rsid"].astype(str))
    pbc = load_pbc(requested_rsids)

    audit_rows = []
    for row in manifest.sort_values("priority").to_dict("records"):
        locus = str(row["pbc_locus"])
        locus_id = f"PBC{int(row['priority']):02d}"
        reference = references[locus]
        merged = reference.merge(pbc, left_on="rsid", right_on="variant_id", how="inner")
        rsid_overlap_n = len(merged)
        merged = merged.loc[
            merged["position_grch37_1based"].eq(merged["base_pair_location"])
        ].copy()
        coordinate_matched_n = len(merged)
        alignments = merged.apply(
            lambda value: align_to_alt(
                value["effect_allele"],
                value["other_allele"],
                value["reference_allele"],
                value["alternate_allele"],
            ),
            axis=1,
        )
        merged["pbc_sign"] = [value[0] for value in alignments]
        merged["pbc_alignment"] = [value[1] for value in alignments]
        alignment_counts = merged["pbc_alignment"].value_counts().to_dict()
        merged = merged.loc[merged["pbc_sign"].notna()].copy()
        merged["pbc_beta_alt"] = merged["beta"] * merged["pbc_sign"]
        merged["pbc_z_alt"] = merged["pbc_beta_alt"] / merged["standard_error"]
        merged = merged.sort_values(["position_grch37_1based", "rsid"]).reset_index(drop=True)
        if len(merged) < 500:
            raise RuntimeError(f"Too few PBC variants for {locus}: {len(merged)}")
        output = merged[
            [
                "rsid",
                "chromosome_grch37",
                "position_grch37_1based",
                "reference_allele",
                "alternate_allele",
                "reference_alt_freq",
                "effect_allele",
                "other_allele",
                "pbc_alignment",
                "pbc_beta_alt",
                "standard_error",
                "p_value",
                "pbc_z_alt",
            ]
        ].rename(
            columns={
                "effect_allele": "pbc_effect_allele_original",
                "other_allele": "pbc_other_allele_original",
                "standard_error": "pbc_se",
                "p_value": "pbc_pvalue",
            }
        )
        output.insert(0, "pbc_locus", locus)
        output.insert(0, "locus_id", locus_id)
        data_path = OUT_DIR / f"{locus_id}__{locus.replace(' ', '_')}.tsv.gz"
        snplist_path = OUT_DIR / f"{locus_id}.snplist"
        output.to_csv(data_path, sep="\t", index=False, compression="gzip")
        snplist_path.write_text("\n".join(output["rsid"]) + "\n", encoding="utf-8")
        audit_rows.append(
            {
                "locus_id": locus_id,
                "priority": row["priority"],
                "pbc_locus": locus,
                "lead_rsid": row["lead_rsid"],
                "reference_panel_variant_n": len(reference),
                "pbc_reference_rsid_overlap_n": rsid_overlap_n,
                "coordinate_matched_n": coordinate_matched_n,
                "harmonized_n": len(output),
                "pbc_alignment_counts": json.dumps(alignment_counts, sort_keys=True),
                "data_path": str(data_path),
                "snplist_path": str(snplist_path),
                "reference_prefix": str(prefixes[locus]),
                "status": "ready",
            }
        )
        print(f"PBC_HARMONIZED\t{locus_id}\t{locus}\t{len(output)}", flush=True)
    pd.DataFrame(audit_rows).to_csv(AUDIT_OUT, sep="\t", index=False)
    print(f"COMPLETE\t{len(audit_rows)}\t{AUDIT_OUT}")


if __name__ == "__main__":
    main()
