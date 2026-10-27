from __future__ import annotations

import csv
import gzip
import math
import re
from pathlib import Path
from statistics import NormalDist


REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data" / "raw" / "genetics" / "finngen" / "r12" / "chirbil_prim" / "regions"
SIGNAL_INPUT = REPO / "work" / "genetics" / "signal_level" / "input"
RESULTS = REPO / "work" / "genetics" / "external_signal_coloc" / "results"

N_CASE = 760
N_CONTROL = 372_273
PRIMARY_ALPHA = 0.05 / 3

LOCUS_CONFIG = {
    "1q23.1": {
        "file": DATA / "01_1q23.1_FCRL3_FCRL3-AS1.tsv.gz",
        "map": SIGNAL_INPUT / "SA05__FCRL3__B_naive.tsv.gz",
        "targets": ["rs945635", "rs3761959"],
        "primary": "rs945635",
    },
    "1p31.3": {
        "file": DATA / "02_1p31.3_IL12RB2.tsv.gz",
        "map": SIGNAL_INPUT / "SA09__IL12RB2__NK.tsv.gz",
        "targets": ["rs6679356"],
        "primary": "rs6679356",
    },
    "18q22.2": {
        "file": DATA / "03_18q22.2_CD226_DOK6.tsv.gz",
        "map": SIGNAL_INPUT / "SA13__CD226__CD4_TCM.tsv.gz",
        "targets": ["rs1808094", "rs17207042", "rs12969657", "rs1790932"],
        "primary": "rs1808094",
    },
}


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rsid_set(value: str) -> set[str]:
    return {item for item in re.split(r"[;,]", value or "") if item and item != "."}


def power_from_expected_z(expected_z: float, alpha: float) -> float:
    normal = NormalDist()
    threshold = normal.inv_cdf(1 - alpha / 2)
    return 1 - normal.cdf(threshold - expected_z) + normal.cdf(-threshold - expected_z)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    variant_rows = []
    locus_rows = []

    for locus, config in LOCUS_CONFIG.items():
        region = read_gzip_tsv(config["file"])
        mapping_rows = read_gzip_tsv(config["map"])
        mapping = {row["rsid"]: row for row in mapping_rows if row["rsid"] in config["targets"]}
        region_pvalues = [float(row["pval"]) for row in region if row["pval"] not in {"", "NA"}]
        top = min(region, key=lambda row: float(row["pval"]))
        locus_rows.append(
            {
                "locus": locus,
                "genome_build": "GRCh38",
                "coordinate_convention": "1-based closed",
                "window_definition": "discovery lead +/- 1 Mb",
                "variant_n": len(region),
                "region_bonferroni_alpha": 0.05 / len(region),
                "top_chromosome": top["#chrom"],
                "top_position_grch38_1based": top["pos"],
                "top_ref": top["ref"],
                "top_alt": top["alt"],
                "top_rsids": top["rsids"],
                "top_beta_alt": top["beta"],
                "top_se": top["sebeta"],
                "top_pvalue": top["pval"],
                "region_bonferroni_significant": float(top["pval"]) < 0.05 / len(region),
                "genome_wide_significant": float(top["pval"]) < 5e-8,
                "file": config["file"].relative_to(REPO).as_posix(),
            }
        )

        for target in config["targets"]:
            discovery = mapping[target]
            expected_position = int(discovery["eqtl_position_grch38_1based"])
            expected_ref = discovery["eqtl_ref_original"]
            expected_alt = discovery["eqtl_alt_original"]
            matches = [row for row in region if target in rsid_set(row["rsids"])]
            if not matches:
                matches = [
                    row
                    for row in region
                    if int(row["pos"]) == expected_position
                    and row["ref"] == expected_ref
                    and row["alt"] == expected_alt
                ]
            if len(matches) != 1:
                variant_rows.append(
                    {
                        "locus": locus,
                        "rsid": target,
                        "is_primary_lead": target == config["primary"],
                        "status": "not_found" if not matches else "ambiguous",
                    }
                )
                continue

            finn = matches[0]
            exact = (
                int(finn["pos"]) == expected_position
                and finn["ref"] == expected_ref
                and finn["alt"] == expected_alt
            )
            if not exact:
                raise ValueError(
                    f"Strict GRCh38 REF/ALT mismatch for {target}: "
                    f"expected {expected_position}:{expected_ref}:{expected_alt}, "
                    f"observed {finn['pos']}:{finn['ref']}:{finn['alt']}"
                )
            pbc_beta = float(discovery["pbc_beta_alt"])
            finn_beta = float(finn["beta"])
            finn_se = float(finn["sebeta"])
            finn_p = float(finn["pval"])
            same_direction = math.copysign(1, pbc_beta) == math.copysign(1, finn_beta)
            one_sided_p = finn_p / 2 if same_direction else 1 - finn_p / 2
            expected_z = abs(pbc_beta) / finn_se
            variant_rows.append(
                {
                    "locus": locus,
                    "rsid": target,
                    "is_primary_lead": target == config["primary"],
                    "status": "exact_match",
                    "genome_build": "GRCh38",
                    "position_grch38_1based": finn["pos"],
                    "ref": finn["ref"],
                    "alt_effect_allele": finn["alt"],
                    "pbc_beta_alt": pbc_beta,
                    "pbc_pvalue": discovery["pbc_pvalue"],
                    "finngen_beta_alt": finn_beta,
                    "finngen_se": finn_se,
                    "finngen_pvalue_two_sided": finn_p,
                    "finngen_pvalue_one_sided_discovery_direction": one_sided_p,
                    "finngen_af_alt": finn["af_alt"],
                    "finngen_af_alt_cases": finn["af_alt_cases"],
                    "finngen_af_alt_controls": finn["af_alt_controls"],
                    "effect_direction_concordant": same_direction,
                    "nominal_replication": same_direction and finn_p < 0.05,
                    "primary_lead_bonferroni_replication": (
                        target == config["primary"] and same_direction and finn_p < PRIMARY_ALPHA
                    ),
                    "expected_z_from_discovery_beta_and_finngen_se": expected_z,
                    "approx_power_two_sided_alpha_0_0167": power_from_expected_z(expected_z, PRIMARY_ALPHA),
                    "n_cases": N_CASE,
                    "n_controls": N_CONTROL,
                }
            )

    write_tsv(
        RESULTS / "04_finngen_region_qc_top_signal.tsv",
        locus_rows,
        list(locus_rows[0].keys()),
    )
    write_tsv(
        RESULTS / "05_finngen_prespecified_variant_replication.tsv",
        variant_rows,
        list(variant_rows[0].keys()),
    )

    primary = [row for row in variant_rows if row.get("is_primary_lead") is True]
    successful = [row for row in primary if row.get("primary_lead_bonferroni_replication") is True]
    nominal = [row for row in primary if row.get("nominal_replication") is True]
    report = f"""# FinnGen R12 PBC 独立 GWAS 复制

## 数据与判定

- 表型：`CHIRBIL_PRIM`（ICD-10 K74.3），760 cases、372,273 controls。
- 坐标：GRCh38，1-based；FinnGen ALT 为效应等位基因。
- 区域：三个发现位点各以发现 lead 的 GRCh38 坐标为中心取 ±1 Mb。
- 主复制位点：rs945635、rs6679356、rs1808094；要求 REF/ALT 精确一致、效应方向一致。
- 名义复制阈值：双侧 P < 0.05；主位点 Bonferroni 阈值：0.05/3 = {PRIMARY_ALPHA:.6g}。

## 结果概览

- 三个位点共提取 {sum(int(row['variant_n']) for row in locus_rows):,} 条区域关联记录。
- 主位点名义复制：{len(nominal)}/3。
- 主位点 Bonferroni 复制：{len(successful)}/3。
- FinnGen 病例数仅 760，未复制不能解释为遗传效应不存在；功效列仅为使用发现效应和 FinnGen 标准误的近似值，可能受赢家诅咒影响。

## 解释边界

- 这是独立 GWAS 复制，不等于外部 eQTL 共定位。
- 多个基因位于同一遗传位点时，只计一个位点级遗传复制证据。
- 公共元数据支持 FinnGen 与 Cordell 2021 主队列高度可能独立，但不能证明逐样本零重叠。
"""
    (RESULTS / "06_finngen_independent_replication_report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
