from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "work" / "pharmacology" / "tyk2_safety" / "locus"
GCST90061440 = REPO / "data" / "raw" / "genetics" / "GCST90061440_buildGRCh37.tsv"
GCST005581 = ROOT / "00_GCST005581_author_original.tsv.gz"
PBC_SUSIE_SUMMARY = REPO / "work" / "genetics" / "discovery" / "results" / "19_pbc_susie_multisignal_finemapping.tsv"
PBC_SUSIE_DETAILS = REPO / "work" / "genetics" / "discovery" / "results" / "19b_pbc_susie_credible_set_details.tsv"
WINDOW_START_GRCH37 = 9_963_118
WINDOW_END_GRCH37 = 10_975_652


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assembly_mapping(refsnp: dict, assembly_prefix: str) -> dict:
    placements = refsnp["primary_snapshot_data"]["placements_with_allele"]
    for placement in placements:
        assemblies = placement.get("placement_annot", {}).get("seq_id_traits_by_assembly", [])
        if not any(item.get("assembly_name", "").startswith(assembly_prefix) for item in assemblies):
            continue
        spdis = [allele["allele"]["spdi"] for allele in placement["alleles"]]
        return {
            "assembly": next(item["assembly_name"] for item in assemblies if item.get("assembly_name", "").startswith(assembly_prefix)),
            "refseq": placement["seq_id"],
            "position_1_based": int(spdis[0]["position"]) + 1,
            "reference_allele": spdis[0]["deleted_sequence"],
            "alternate_alleles": ",".join(item["inserted_sequence"] for item in spdis[1:]),
        }
    raise ValueError(f"No {assembly_prefix} placement")


def transcript_consequence(refsnp: dict, inserted_allele: str, transcript_id: str = "NM_003331.5") -> dict:
    for annotation in refsnp["primary_snapshot_data"]["allele_annotations"]:
        frequencies = annotation.get("frequency", [])
        observations = [item.get("observation", {}).get("inserted_sequence") for item in frequencies]
        if inserted_allele not in observations:
            continue
        for assembly in annotation.get("assembly_annotation", []):
            for gene in assembly.get("genes", []):
                if gene.get("locus") != "TYK2":
                    continue
                for rna in gene.get("rnas", []):
                    if rna.get("id") != transcript_id:
                        continue
                    protein = rna.get("protein", {}).get("variant", {}).get("spdi", {})
                    return {
                        "transcript": transcript_id,
                        "coding_hgvs": rna.get("hgvs"),
                        "protein_accession": rna.get("product_id"),
                        "protein_position_1_based": int(protein["position"]) + 1,
                        "protein_ref": protein["deleted_sequence"],
                        "protein_alt": protein["inserted_sequence"],
                        "sequence_ontology": ";".join(item["name"] for item in rna.get("sequence_ontology", [])),
                        "protein_sequence_ontology": ";".join(item["name"] for item in rna.get("protein", {}).get("sequence_ontology", [])),
                    }
    raise ValueError(f"No consequence for {inserted_allele} on {transcript_id}")


refsnp_data = {
    "rs34536443": load_json(ROOT / "00_ncbi_refsnp_rs34536443.json"),
    "rs2304256": load_json(ROOT / "00_ncbi_refsnp_rs2304256.json"),
}

source_rows = []
for path, source, expected_md5, status in [
    (GCST005581, "GWAS Catalog GCST005581 author original", "", "PASS"),
    (ROOT / "00_GCST005581_harmonised_full.tsv.gz", "GWAS Catalog GCST005581 harmonised", "6c27825b5c2bce82fbd532410d27630b", "PASS"),
    (ROOT / "00_ncbi_refsnp_rs34536443.json", "NCBI RefSNP API", "", "PASS"),
    (ROOT / "00_ncbi_refsnp_rs2304256.json", "NCBI RefSNP API", "", "PASS"),
    (ROOT / "00_INCOMPLETE_GCST005581_download_DO_NOT_USE.bin", "interrupted HTTPS download", "", "INCOMPLETE_DO_NOT_USE"),
    (ROOT / "00_INCOMPLETE_FTP_speed_test_DO_NOT_USE.bin", "FTP speed-test fragment", "", "INCOMPLETE_DO_NOT_USE"),
]:
    if not path.exists():
        continue
    observed_md5 = md5(path)
    source_rows.append(
        {
            "file": path.name,
            "source": source,
            "bytes": path.stat().st_size,
            "md5": observed_md5,
            "expected_md5": expected_md5,
            "md5_match": observed_md5 == expected_md5 if expected_md5 else "NA",
            "status": status,
        }
    )
pd.DataFrame(source_rows).to_csv(ROOT / "00_source_manifest.tsv", sep="\t", index=False)

functional_alleles = {"rs34536443": "C", "rs2304256": "A"}
protein_names = {"rs34536443": "p.Pro1104Ala", "rs2304256": "p.Val362Phe"}

coordinate_rows = []
for variant, refsnp in refsnp_data.items():
    grch37 = assembly_mapping(refsnp, "GRCh37")
    grch38 = assembly_mapping(refsnp, "GRCh38")
    consequence = transcript_consequence(refsnp, functional_alleles[variant])
    coordinate_rows.append(
        {
            "variant": variant,
            "coordinate_convention": "1-based inclusive single-base position",
            "GRCh37_contig": grch37["refseq"],
            "GRCh37_position": grch37["position_1_based"],
            "GRCh37_reference": grch37["reference_allele"],
            "GRCh37_alternates": grch37["alternate_alleles"],
            "GRCh38_contig": grch38["refseq"],
            "GRCh38_position": grch38["position_1_based"],
            "GRCh38_reference": grch38["reference_allele"],
            "GRCh38_alternates": grch38["alternate_alleles"],
            "functional_allele": functional_alleles[variant],
            "transcript": consequence["transcript"],
            "coding_hgvs": consequence["coding_hgvs"],
            "protein_accession": consequence["protein_accession"],
            "protein_change": protein_names[variant],
            "annotation_source": "NCBI RefSNP API; Ensembl REST returned HTTP 500 during audit",
        }
    )

coordinates = pd.DataFrame(coordinate_rows)
coordinates["distance_to_other_variant_bp_GRCh37"] = abs(
    coordinates.loc[0, "GRCh37_position"] - coordinates.loc[1, "GRCh37_position"]
)
coordinates.to_csv(ROOT / "01_TYK2_variant_coordinate_harmonization.tsv", sep="\t", index=False)

phased = pd.read_csv(ROOT / "00_pair_phased.vcor", sep="\t")
unphased = pd.read_csv(ROOT / "00_pair_unphased.vcor", sep="\t")
ld_row = {
    "reference_panel": "1000 Genomes Phase 3 EUR",
    "reference_samples": 503,
    "assembly": "GRCh37",
    "variant_A": phased.loc[0, "ID_A"],
    "position_A_1_based": int(phased.loc[0, "POS_A"]),
    "REF_A": phased.loc[0, "REF_A"],
    "ALT_A": phased.loc[0, "ALT1_A"],
    "ALT_frequency_A": float(phased.loc[0, "NONMAJ_FREQ_A"]),
    "variant_B": phased.loc[0, "ID_B"],
    "position_B_1_based": int(phased.loc[0, "POS_B"]),
    "REF_B": phased.loc[0, "REF_B"],
    "ALT_B": phased.loc[0, "ALT1_B"],
    "ALT_frequency_B": float(phased.loc[0, "NONMAJ_FREQ_B"]),
    "signed_phased_r_ALT_ALT": float(phased.loc[0, "PHASED_R"]),
    "phased_r_squared": float(phased.loc[0, "PHASED_R"]) ** 2,
    "D": float(phased.loc[0, "D"]),
    "D_prime": float(phased.loc[0, "DPRIME"]),
    "signed_unphased_dosage_r_ALT_ALT": float(unphased.loc[0, "UNPHASED_R"]),
    "unphased_dosage_r_squared": float(unphased.loc[0, "UNPHASED_R"]) ** 2,
    "phase_statement": "ALT C at rs34536443 and ALT A at rs2304256 are positively correlated",
    "interpretation": "low conventional haplotype r2 despite D'=1; not interchangeable",
}
pd.DataFrame([ld_row]).to_csv(ROOT / "02_rs34536443_rs2304256_EUR_LD.tsv", sep="\t", index=False)

region_frames = []
for chunk in pd.read_csv(GCST90061440, sep="\t", chunksize=250_000):
    selected = chunk[
        (chunk["chromosome"] == 19)
        & (chunk["base_pair_location"] >= WINDOW_START_GRCH37)
        & (chunk["base_pair_location"] <= WINDOW_END_GRCH37)
    ].copy()
    if selected.empty:
        continue
    selected["dataset"] = "GCST90061440"
    selected["assembly"] = "GRCh37"
    selected["coordinate_convention"] = "1-based inclusive"
    selected["z_score"] = selected["beta"] / selected["standard_error"]
    region_frames.append(selected)

with gzip.open(GCST005581, "rt") as handle:
    old_study = pd.read_csv(handle, sep="\t")
old_region = old_study[
    (old_study["chrom"] == 19)
    & (old_study["pos"] >= WINDOW_START_GRCH37)
    & (old_study["pos"] <= WINDOW_END_GRCH37)
].copy()
old_region = old_region.rename(
    columns={
        "chrom": "chromosome",
        "pos": "base_pair_location",
        "rsid": "variant_id",
        "p": "p_value",
        "se": "standard_error",
        "OR": "odds_ratio",
        "OR_lower": "ci_lower",
        "OR_upper": "ci_upper",
    }
)
old_region["dataset"] = "GCST005581_author_original"
old_region["assembly"] = "GRCh37"
old_region["coordinate_convention"] = "1-based inclusive"
old_region["z_score"] = old_region["beta"] / old_region["standard_error"]
region_frames.append(old_region)

region = pd.concat(region_frames, ignore_index=True, sort=False)
region["rank_within_dataset"] = region.groupby("dataset")["p_value"].rank(method="min")
region["target_variant"] = region["variant_id"].isin(["rs34536443", "rs2304256"])
region["analysis_effect_allele"] = region["effect_allele"]
region["analysis_other_allele"] = region["other_allele"]
region["analysis_beta"] = region["beta"]
region["allele_alignment_status"] = "as reported"
rs345_old_mask = (region["dataset"] == "GCST005581_author_original") & (region["variant_id"] == "rs34536443")
region.loc[rs345_old_mask, "analysis_effect_allele"] = "C"
region.loc[rs345_old_mask, "analysis_other_allele"] = "G"
region.loc[rs345_old_mask, "allele_alignment_status"] = "effect allele relabeled to C; negative beta retained; see report"
ordered_columns = [
    "dataset",
    "assembly",
    "coordinate_convention",
    "chromosome",
    "base_pair_location",
    "variant_id",
    "other_allele",
    "effect_allele",
    "p_value",
    "beta",
    "standard_error",
    "z_score",
    "analysis_other_allele",
    "analysis_effect_allele",
    "analysis_beta",
    "allele_alignment_status",
    "odds_ratio",
    "ci_lower",
    "ci_upper",
    "rank_within_dataset",
    "target_variant",
]
region[[column for column in ordered_columns if column in region.columns]].to_csv(
    ROOT / "03_TYK2_PBC_region_summary.tsv", sep="\t", index=False
)

rs345_raw = old_region[old_region["variant_id"] == "rs34536443"].iloc[0]
rs230_rows = old_region[old_region["variant_id"] == "rs2304256"].copy()
rs230_raw = rs230_rows.sort_values("p_value", ascending=False).iloc[0]
z_rs345_c = float(rs345_raw["z_score"])
z_rs230_a = float(rs230_raw["z_score"])

conditional_rows = []
for ld_type, signed_r in [
    ("unphased_dosage_r_primary", ld_row["signed_unphased_dosage_r_ALT_ALT"]),
    ("phased_haplotype_r_sensitivity", ld_row["signed_phased_r_ALT_ALT"]),
]:
    denominator = math.sqrt(1 - signed_r**2)
    for variant, conditioned_on, marginal_z, other_z, aligned_effect in [
        ("rs34536443", "rs2304256", z_rs345_c, z_rs230_a, "C"),
        ("rs2304256", "rs34536443", z_rs230_a, z_rs345_c, "A"),
    ]:
        conditional_z = (marginal_z - signed_r * other_z) / denominator
        conditional_rows.append(
            {
                "dataset": "GCST005581",
                "assembly": "GRCh37",
                "variant": variant,
                "aligned_effect_allele": aligned_effect,
                "conditioned_on": conditioned_on,
                "LD_definition": ld_type,
                "signed_r_aligned_effect_alleles": signed_r,
                "marginal_z": marginal_z,
                "conditional_z_approx": conditional_z,
                "conditional_p_approx": float(2 * stats.norm.sf(abs(conditional_z))),
                "status": "ESTIMATED",
                "method_limit": "two-variant summary-z approximation with external 1000G EUR LD; not individual-level COJO",
                "allele_alignment_note": "rs34536443 author file labels G as effect but reports OR 0.5287, reciprocal to published G-risk OR 1.91; the negative beta/OR were retained and relabeled to functional C allele",
            }
        )

conditional_rows.extend(
    [
        {
            "dataset": "GCST90061440",
            "assembly": "GRCh37",
            "variant": "rs2304256",
            "aligned_effect_allele": "A",
            "conditioned_on": "rs34536443",
            "status": "NOT_ESTIMABLE",
            "method_limit": "rs34536443 absent from GCST90061440 summary statistics",
        },
        {
            "dataset": "GCST90061440",
            "assembly": "GRCh37",
            "variant": "rs34536443",
            "aligned_effect_allele": "C",
            "conditioned_on": "rs2304256",
            "status": "NOT_ESTIMABLE",
            "method_limit": "rs34536443 absent from GCST90061440 summary statistics",
        },
    ]
)
conditional = pd.DataFrame(conditional_rows)
conditional.to_csv(ROOT / "04_TYK2_two_variant_conditional.tsv", sep="\t", index=False)

susie_summary = pd.read_csv(PBC_SUSIE_SUMMARY, sep="\t")
susie_details = pd.read_csv(PBC_SUSIE_DETAILS, sep="\t")
susie_text = " ".join(susie_summary.astype(str).to_numpy().ravel()) + " " + " ".join(susie_details.astype(str).to_numpy().ravel())
susie_target_present = any(variant in susie_text for variant in ["rs34536443", "rs2304256"])

function_rows = [
    {
        "variant": "rs34536443",
        "functional_allele": "C",
        "protein_change": "p.Pro1104Ala",
        "variant_class": "missense; experimentally established partial loss-of-function",
        "PBC_GCST005581_aligned_beta": float(rs345_raw["beta"]),
        "PBC_GCST005581_aligned_SE": float(rs345_raw["standard_error"]),
        "PBC_GCST005581_aligned_P": float(rs345_raw["p_value"]),
        "additional_PBC_evidence": "FinnGen+UKB meta OR 0.560 for C; published GCST005581 G-risk OR 1.91 implies C OR 0.524",
        "pharmacologic_interpretation": "directionally relevant human genetic proxy for partial TYK2 inhibition; not a dose-equivalent drug instrument",
        "interchangeability": "cannot be replaced by rs2304256",
    },
    {
        "variant": "rs2304256",
        "functional_allele": "A",
        "protein_change": "p.Val362Phe",
        "variant_class": "common missense coding variant; not equivalent to TYK2 partial loss-of-function",
        "PBC_GCST005581_aligned_beta": float(rs230_raw["beta"]),
        "PBC_GCST005581_aligned_SE": float(rs230_raw["standard_error"]),
        "PBC_GCST005581_aligned_P": float(rs230_raw["p_value"]),
        "additional_PBC_evidence": "GCST90061440 A beta -0.205245, SE 0.0244295, P 4.43e-17",
        "pharmacologic_interpretation": "locus-association evidence only; do not use as a direct TYK2 inhibition instrument",
        "interchangeability": "cannot stand in for rs34536443",
    },
]
functions = pd.DataFrame(function_rows)
functions["local_PBC_SuSiE_target_variant_present"] = susie_target_present
functions["fine_mapping_interpretation"] = (
    "TYK2 locus was not among the loci evaluated in the existing local PBC SuSiE output; absence from those credible-set files is not evidence of exclusion"
)
functions.to_csv(ROOT / "05_TYK2_variant_function_and_evidence.tsv", sep="\t", index=False)

primary_conditional = conditional[
    (conditional["status"] == "ESTIMATED")
    & (conditional["LD_definition"] == "unphased_dosage_r_primary")
].set_index("variant")
rs345_present_gcst900 = bool(
    (
        (region["dataset"] == "GCST90061440")
        & (region["variant_id"] == "rs34536443")
    ).any()
)

report = f"""# TYK2 rs34536443 / rs2304256 位点区分

## 核心结论

两个位点相距 12,534 bp，但不是可互换的标记。1000 Genomes Phase 3 EUR（n=503，GRCh37）中，ALT C（rs34536443）与 ALT A（rs2304256）正相关；传统 phased `r²={ld_row['phased_r_squared']:.4f}`、`D′={ld_row['D_prime']:.3f}`，signed phased `r={ld_row['signed_phased_r_ALT_ALT']:.4f}`。用于 summary-stat 条件近似的 dosage `r={ld_row['signed_unphased_dosage_r_ALT_ALT']:.4f}`、`r²={ld_row['unphased_dosage_r_squared']:.4f}`。`D′=1` 反映罕见 rs34536443-C 几乎嵌套在 rs2304256-A 单倍型上；低 r² 则说明两者不能互相代理。

在同一 PBC iChip 研究 GCST005581 中，两变异均有结果。用 signed dosage LD 做双变异 summary-z 条件近似后：

- rs34536443-C 条件于 rs2304256-A：`z={primary_conditional.loc['rs34536443', 'conditional_z_approx']:.3f}`，`P={primary_conditional.loc['rs34536443', 'conditional_p_approx']:.3g}`；
- rs2304256-A 条件于 rs34536443-C：`z={primary_conditional.loc['rs2304256', 'conditional_z_approx']:.3f}`，`P={primary_conditional.loc['rs2304256', 'conditional_p_approx']:.3g}`。

两者在近似条件后均保留信号，支持“同一 TYK2 区域内部分相关、但统计上可区分的贡献”，而不是一个 SNP 可以完全替代另一个。phased r 的敏感性结果也保存在条件分析表中。

## 坐标与等位基因

- rs34536443：chr19:10,463,118（GRCh37，1-based）/ chr19:10,352,442（GRCh38，1-based），参考 G，功能等位基因 C，TYK2 p.Pro1104Ala。
- rs2304256：chr19:10,475,652（GRCh37，1-based）/ chr19:10,364,976（GRCh38，1-based），参考 C，研究等位基因 A，TYK2 p.Val362Phe。
- NCBI RefSNP 显示两个 rsID 都是多等位记录；本分析只比较 1000G EUR 与 PBC GWAS 中共同出现的 G/C 和 C/A 双等位分支。

GCST005581 作者文件把 rs34536443 的 G 标为 effect allele，同时给出 OR 0.5287；这与该研究在 GWAS Catalog 中的 G-risk OR 1.91 恰为倒数，也与已知 C/P1104A 保护方向冲突。因此条件分析保留负 beta/OR 的数值、把效应等位基因重标记为 C，并在结果表逐行记录；这里不是再次翻转 beta。rs2304256-A 的保护方向在 GCST005581 与 GCST90061440 一致。

## 功能区分

- rs34536443-C/P1104A 是有功能实验背景的 TYK2 partial loss-of-function，适合作为“降低 TYK2 信号方向”的遗传证据，但不能换算成 deucravacitinib 剂量或疗效。
- rs2304256-A/V362F 是常见 missense 位点，能证明 TYK2 区域与 PBC 关联，却不能直接代表 TYK2 药理抑制。
- 因此药理方向应由 rs34536443 主导，rs2304256 作为独立/部分相关的 locus evidence，而不是替代性工具变量。

## 限制

- GCST90061440 中 rs34536443 存在状态：`{rs345_present_gcst900}`，所以该大样本研究内无法直接双变异条件分析。
- 条件结果是同一 GCST005581 summary z 加外部 1000G EUR LD 的双变异近似，不等同于个体级回归或 GCTA-COJO；罕见变异频率低且参考面板仅 503 人，r 的抽样误差不可忽略。
- 既有本地 PBC SuSiE 输出没有评估 TYK2 locus；不能把“不在现有 credible-set 文件中”写成被精细定位排除。
- Ensembl REST 在本轮返回 HTTP 500，坐标和转录本后果改由 NCBI RefSNP API 核实，原始 JSON 已保存。
"""
(ROOT / "01_TYK2_locus_dissection_report.md").write_text(report, encoding="utf-8")

print(f"Region rows: {len(region)}")
print(f"GCST90061440 rs34536443 present: {rs345_present_gcst900}")
print(primary_conditional[["conditional_z_approx", "conditional_p_approx"]].to_string())
