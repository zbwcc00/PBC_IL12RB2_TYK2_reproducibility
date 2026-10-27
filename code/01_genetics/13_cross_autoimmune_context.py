from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
WORK_DIR = REPO / "work" / "genetics" / "cross_autoimmune_context"
REGION_DIR = WORK_DIR / "raw" / "regions"
PROCESSED_DIR = WORK_DIR / "processed"
REPORT_DIR = WORK_DIR / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)

PBC_INPUT = REPO / "work" / "genetics" / "fine_mapping" / "input" / "PBC04__1p31.3.tsv.gz"
LD_PREFIX = REPO / "work" / "genetics" / "fine_mapping" / "ld" / "PBC04.unphased.vcor1.bin"
PBC_CS_PATH = (
    REPO
    / "work"
    / "genetics"
    / "regulatory_annotation"
    / "processed"
    / "01_credible_set_variant_regulatory_annotation.tsv"
)

STUDY_META = {
    "PSC": {"accession": "GCST004030", "cases": 2871, "controls": 12019},
    "RA": {"accession": "GCST90837201", "cases": 26251, "controls": 431261},
    "IBD": {"accession": "GCST004131", "cases": 25042, "controls": 34915},
    "Psoriasis": {"accession": "GCST90472771", "cases": 36466, "controls": 458078},
}
COMPONENT_LEADS = {"L1": "rs6679356", "L2": "rs3828069", "L3": "rs112588573"}
COMPLEMENTS = str.maketrans("ACGT", "TGCA")


def write_tsv(path: Path, frame: pd.DataFrame):
    frame.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def logsumexp(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if not len(finite):
        return -math.inf
    maximum = float(np.max(finite))
    return maximum + math.log(float(np.exp(finite - maximum).sum()))


def logdiffexp(left: float, right: float) -> float:
    if right >= left:
        return -math.inf
    return left + math.log1p(-math.exp(right - left))


def log_abf(beta: np.ndarray, variance: np.ndarray, prior_variance: float = 0.04) -> np.ndarray:
    shrinkage = prior_variance / (prior_variance + variance)
    z_squared = beta**2 / variance
    return 0.5 * (np.log1p(-shrinkage) + shrinkage * z_squared)


def coloc_abf(
    log_bf_1: np.ndarray,
    log_bf_2: np.ndarray,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 1e-5,
) -> dict:
    sum_1 = logsumexp(log_bf_1)
    sum_2 = logsumexp(log_bf_2)
    sum_12 = logsumexp(log_bf_1 + log_bf_2)
    sum_distinct = logdiffexp(sum_1 + sum_2, sum_12)
    log_hypotheses = np.array(
        [
            0.0,
            math.log(p1) + sum_1,
            math.log(p2) + sum_2,
            math.log(p1) + math.log(p2) + sum_distinct,
            math.log(p12) + sum_12,
        ]
    )
    denominator = logsumexp(log_hypotheses)
    posterior = np.exp(log_hypotheses - denominator)
    return {f"PP.H{index}": float(value) for index, value in enumerate(posterior)}


def complement(allele: str) -> str:
    if not allele or any(base not in "ACGT" for base in allele):
        return ""
    return allele.translate(COMPLEMENTS)


def align_beta_to_alt(
    effect_allele: str,
    other_allele: str,
    beta: float,
    reference: str,
    alternate: str,
) -> tuple[float, str]:
    effect_allele = effect_allele.upper()
    other_allele = other_allele.upper()
    reference = reference.upper()
    alternate = alternate.upper()
    if effect_allele == alternate and other_allele == reference:
        return beta, "direct_alt"
    if effect_allele == reference and other_allele == alternate:
        return -beta, "flipped_ref"
    effect_complement = complement(effect_allele)
    other_complement = complement(other_allele)
    if effect_complement == alternate and other_complement == reference:
        return beta, "strand_direct_alt"
    if effect_complement == reference and other_complement == alternate:
        return -beta, "strand_flipped_ref"
    return math.nan, "allele_mismatch"


def standardize_trait(path: Path, trait: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "rsid" in frame.columns:
        rsid = frame["rsid"]
    elif "hm_rsid" in frame.columns:
        rsid = frame["hm_rsid"]
    else:
        rsid = frame["variant_id"]
    beta = pd.to_numeric(frame.get("beta", ""), errors="coerce")
    if beta.isna().all() and "odds_ratio" in frame.columns:
        beta = np.log(pd.to_numeric(frame["odds_ratio"], errors="coerce"))
    standardized = pd.DataFrame(
        {
            "trait": trait,
            "rsid": rsid,
            "chromosome": frame["chromosome"].str.removeprefix("chr"),
            "position_grch38": pd.to_numeric(frame["base_pair_location"], errors="coerce"),
            "effect_allele": frame["effect_allele"].str.upper(),
            "other_allele": frame["other_allele"].str.upper(),
            "beta": beta,
            "se": pd.to_numeric(frame["standard_error"], errors="coerce"),
            "pvalue": pd.to_numeric(frame["p_value"], errors="coerce"),
            "eaf": pd.to_numeric(frame.get("effect_allele_frequency", ""), errors="coerce"),
        }
    )
    standardized = standardized[
        standardized["rsid"].str.startswith("rs")
        & standardized["beta"].notna()
        & standardized["se"].gt(0)
        & standardized["pvalue"].notna()
    ].copy()
    standardized = standardized.sort_values("pvalue").drop_duplicates("rsid", keep="first")
    return standardized


def load_pbc() -> pd.DataFrame:
    pbc = pd.read_csv(PBC_INPUT, sep="\t", compression="gzip", dtype=str).fillna("")
    for column in ("pbc_beta_alt", "pbc_se", "pbc_pvalue"):
        pbc[column] = pd.to_numeric(pbc[column], errors="coerce")
    pbc = pbc[pbc["rsid"].str.startswith("rs") & pbc["pbc_se"].gt(0)].copy()
    pbc = pbc.sort_values("pbc_pvalue").drop_duplicates("rsid", keep="first")
    return pbc


def load_ld() -> tuple[list[str], np.ndarray]:
    variants_path = Path(f"{LD_PREFIX}.vars")
    variants = variants_path.read_text(encoding="utf-8").splitlines()
    matrix = np.fromfile(LD_PREFIX, dtype=np.float32)
    expected = len(variants) ** 2
    if matrix.size != expected:
        raise RuntimeError(f"LD size mismatch: {matrix.size} versus {expected}")
    return variants, matrix.reshape((len(variants), len(variants)))


def main():
    pbc = load_pbc()
    pbc_cs = pd.read_csv(PBC_CS_PATH, sep="\t", dtype=str).fillna("")
    pbc_cs_lookup = pbc_cs.set_index("rsid").to_dict("index")
    ld_variants, ld_matrix = load_ld()
    ld_index = {variant: index for index, variant in enumerate(ld_variants)}

    summary_rows = []
    coloc_rows = []
    component_rows = []
    sentinel_rows = []
    harmonization_rows = []

    for trait, study in STUDY_META.items():
        trait_path = REGION_DIR / f"{trait}_IL12RB2_GRCh38.tsv"
        trait_frame = standardize_trait(trait_path, trait)
        lead = trait_frame.nsmallest(1, "pvalue").iloc[0]

        merged = pbc.merge(trait_frame, on="rsid", how="inner", suffixes=("_pbc", "_trait"))
        aligned_values = []
        alignments = []
        for _, row in merged.iterrows():
            aligned, alignment = align_beta_to_alt(
                row["effect_allele"],
                row["other_allele"],
                float(row["beta"]),
                row["reference_allele"],
                row["alternate_allele"],
            )
            aligned_values.append(aligned)
            alignments.append(alignment)
        merged["trait_beta_alt"] = aligned_values
        merged["alignment"] = alignments
        harmonization_rows.append(
            {
                "trait": trait,
                "accession": study["accession"],
                "regional_rows": len(trait_frame),
                "rsid_overlap_with_pbc": len(merged),
                "allele_aligned_rows": int(merged["trait_beta_alt"].notna().sum()),
                "allele_mismatch_rows": int(merged["trait_beta_alt"].isna().sum()),
            }
        )
        merged = merged[merged["trait_beta_alt"].notna()].copy()
        merged = merged[merged["rsid"].isin(ld_index)].copy()
        merged["pbc_variance"] = merged["pbc_se"] ** 2
        merged["trait_variance"] = merged["se"] ** 2
        merged["pbc_log_abf"] = log_abf(
            merged["pbc_beta_alt"].to_numpy(float),
            merged["pbc_variance"].to_numpy(float),
        )
        merged["trait_log_abf"] = log_abf(
            merged["trait_beta_alt"].to_numpy(float),
            merged["trait_variance"].to_numpy(float),
        )

        trait_denominator = logsumexp(merged["trait_log_abf"].to_numpy(float))
        merged["trait_abf_posterior"] = np.exp(merged["trait_log_abf"] - trait_denominator)
        pbc_denominator = logsumexp(merged["pbc_log_abf"].to_numpy(float))
        merged["pbc_abf_posterior"] = np.exp(merged["pbc_log_abf"] - pbc_denominator)

        coloc_default = None
        for p12 in (1e-6, 1e-5, 1e-4):
            posterior = coloc_abf(
                merged["pbc_log_abf"].to_numpy(float),
                merged["trait_log_abf"].to_numpy(float),
                p12=p12,
            )
            coloc_rows.append(
                {
                    "trait": trait,
                    "accession": study["accession"],
                    "overlap_variants": len(merged),
                    "p1": 1e-4,
                    "p2": 1e-4,
                    "p12": p12,
                    **posterior,
                }
            )
            if p12 == 1e-5:
                coloc_default = posterior

        lead_r2 = {}
        if lead["rsid"] in ld_index:
            lead_index = ld_index[lead["rsid"]]
            for component, component_lead in COMPONENT_LEADS.items():
                lead_r2[component] = float(ld_matrix[lead_index, ld_index[component_lead]] ** 2)
        else:
            lead_r2 = {component: math.nan for component in COMPONENT_LEADS}

        for component, component_lead in COMPONENT_LEADS.items():
            component_index = ld_index[component_lead]
            variant_indices = np.array([ld_index[variant] for variant in merged["rsid"]])
            r_squared = ld_matrix[variant_indices, component_index].astype(float) ** 2
            component_mass = float(merged.loc[r_squared >= 0.8, "trait_abf_posterior"].sum())
            component_rows.append(
                {
                    "trait": trait,
                    "accession": study["accession"],
                    "component": component,
                    "component_lead": component_lead,
                    "trait_lead_rsid": lead["rsid"],
                    "trait_lead_pvalue": lead["pvalue"],
                    "trait_lead_r2_to_component": lead_r2[component],
                    "trait_abf_posterior_mass_r2_ge_0_8": component_mass,
                    "variants_in_component_proxy_set": int((r_squared >= 0.8).sum()),
                }
            )

        region_significant = bool(lead["pvalue"] < 5e-8)
        max_component_r2 = max(value for value in lead_r2.values() if math.isfinite(value)) if any(
            math.isfinite(value) for value in lead_r2.values()
        ) else math.nan
        max_component_mass = max(
            row["trait_abf_posterior_mass_r2_ge_0_8"]
            for row in component_rows
            if row["trait"] == trait
        )
        if region_significant and coloc_default["PP.H4"] >= 0.5 and max_component_r2 >= 0.8:
            classification = "shared_causal_signal"
        elif region_significant and coloc_default["PP.H3"] >= 0.8 and max_component_r2 < 0.2:
            classification = "distinct_local_component"
        elif not region_significant:
            classification = "no_detectable_regional_association"
        else:
            classification = "indeterminate"

        summary_rows.append(
            {
                "trait": trait,
                "accession": study["accession"],
                "cases": study["cases"],
                "controls": study["controls"],
                "regional_variants": len(trait_frame),
                "overlap_variants_for_coloc": len(merged),
                "regional_lead_rsid": lead["rsid"],
                "regional_lead_position_grch38": int(lead["position_grch38"]),
                "regional_lead_beta": lead["beta"],
                "regional_lead_se": lead["se"],
                "regional_lead_pvalue": lead["pvalue"],
                "genome_wide_significant_region": region_significant,
                "lead_r2_to_PBC_L1": lead_r2["L1"],
                "lead_r2_to_PBC_L2": lead_r2["L2"],
                "lead_r2_to_PBC_L3": lead_r2["L3"],
                "max_PBC_component_posterior_mass": max_component_mass,
                "coloc_PP_H3_default": coloc_default["PP.H3"],
                "coloc_PP_H4_default": coloc_default["PP.H4"],
                "classification": classification,
            }
        )

        for sentinel in pbc_cs["rsid"]:
            sentinel_match = merged[merged["rsid"] == sentinel]
            if sentinel_match.empty:
                sentinel_rows.append(
                    {
                        "trait": trait,
                        "accession": study["accession"],
                        "rsid": sentinel,
                        "pbc_component": pbc_cs_lookup[sentinel]["component"],
                        "available": False,
                    }
                )
                continue
            row = sentinel_match.iloc[0]
            pbc_effect_allele = pbc_cs_lookup[sentinel]["gwas_effect_allele"]
            pbc_other_allele = pbc_cs_lookup[sentinel]["gwas_other_allele"]
            pbc_beta_effect = float(pbc_cs_lookup[sentinel]["gwas_beta"])
            if pbc_effect_allele == row["alternate_allele"]:
                trait_beta_effect = float(row["trait_beta_alt"])
            elif pbc_effect_allele == row["reference_allele"]:
                trait_beta_effect = -float(row["trait_beta_alt"])
            else:
                trait_beta_effect = math.nan
            risk_allele = pbc_effect_allele if pbc_beta_effect > 0 else pbc_other_allele
            trait_beta_risk = trait_beta_effect if pbc_beta_effect > 0 else -trait_beta_effect
            sentinel_rows.append(
                {
                    "trait": trait,
                    "accession": study["accession"],
                    "rsid": sentinel,
                    "pbc_component": pbc_cs_lookup[sentinel]["component"],
                    "available": True,
                    "pbc_effect_allele": pbc_effect_allele,
                    "pbc_log_odds_per_effect_allele": pbc_beta_effect,
                    "pbc_risk_increasing_allele": risk_allele,
                    "trait_beta_per_pbc_effect_allele": trait_beta_effect,
                    "trait_beta_per_pbc_risk_allele": trait_beta_risk,
                    "trait_se": row["se"],
                    "trait_pvalue": row["pvalue"],
                    "trait_95ci_lower_per_pbc_risk_allele": trait_beta_risk - 1.96 * float(row["se"]),
                    "trait_95ci_upper_per_pbc_risk_allele": trait_beta_risk + 1.96 * float(row["se"]),
                }
            )

        write_tsv(
            PROCESSED_DIR / f"04_{trait}_PBC_harmonized_coloc_input.tsv",
            merged,
        )

    summary = pd.DataFrame(summary_rows)
    write_tsv(PROCESSED_DIR / "05_cross_trait_signal_context_summary.tsv", summary)
    write_tsv(PROCESSED_DIR / "06_cross_trait_coloc_abf_sensitivity.tsv", pd.DataFrame(coloc_rows))
    write_tsv(PROCESSED_DIR / "07_cross_trait_component_posterior_mass.tsv", pd.DataFrame(component_rows))
    write_tsv(PROCESSED_DIR / "08_cross_trait_PBC_sentinel_effects.tsv", pd.DataFrame(sentinel_rows))
    write_tsv(PROCESSED_DIR / "09_cross_trait_harmonization_qc.tsv", pd.DataFrame(harmonization_rows))

    report_lines = [
        "# Step 3：IL12RB2 区域跨自身免疫病 signal-level context map",
        "",
        "## 判定",
        "",
        "**GO（作为疾病背景定位图），但结论是‘同一区域存在不同免疫病信号’，而不是泛自身免疫共享。**",
        "",
        "四个疾病均使用公开、欧洲祖源优先、GRCh38协调化区域统计。分类同时依据区域显著性、1000 Genomes EUR LD、Wakefield ABF 共定位以及 PBC L1/L2/L3 proxy 集合内的后验质量。",
        "",
        "## 结果",
        "",
    ]
    for _, row in summary.iterrows():
        report_lines.append(
            f"- **{row['trait']} ({row['accession']})**：区域 lead `{row['regional_lead_rsid']}`，"
            f"P={float(row['regional_lead_pvalue']):.3g}；与 PBC L1/L2/L3 的 r²分别为 "
            f"{float(row['lead_r2_to_PBC_L1']):.3f}/{float(row['lead_r2_to_PBC_L2']):.3f}/{float(row['lead_r2_to_PBC_L3']):.3f}；"
            f"PP.H3={float(row['coloc_PP_H3_default']):.3f}，PP.H4={float(row['coloc_PP_H4_default']):.3f}；"
            f"分类：`{row['classification']}`。"
        )
    report_lines.extend(
        [
            "",
            "## 生物学解释",
            "",
            "IBD 与 psoriasis 的强峰由彼此高度相关的 IL23R 邻近 leads 驱动，但这些 leads 与 PBC 的三个精细定位成分均处于低 LD，且全区域共定位支持 H3 而非 H4。这说明 1p31.3 是一个多疾病免疫调控邻域，但 PBC 的 IL12RB2 成分并不是 IBD/psoriasis IL23R 峰的简单替代标记。",
            "",
            "RA 与 PSC 在本区域未达到全基因组显著，故不能用‘阴性’证明 PBC 特异性。更准确的结论是：在当前样本量和表型定义下，没有检测到与 PBC L1/L2 相容的区域信号。PSC 数据的病例数明显较少，应保留功效限制；RA 样本量较大，其结果对排除同等量级效应更有信息性。",
            "",
            "## 稿件价值",
            "",
            "该 context map 提升了主线的定位精度：PBC L1 继续承担主要遗传锚定，PBC L2 连接激活记忆 CD4 T-cell IL12RB2 调控，而 IBD/psoriasis 则定义了邻近但独立的 IL23R 自身免疫成分。结果避免把整个 1p31.3 区域笼统解释为共享 IL12RB2 信号。",
        ]
    )
    (REPORT_DIR / "03_step3_cross_autoimmune_context_report_CN.md").write_text(
        "\n".join(report_lines) + "\n",
        encoding="utf-8",
    )

    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": {
            "coloc": "Wakefield ABF single-causal-variant screen with p1=p2=1e-4 and p12 sensitivity 1e-6, 1e-5, 1e-4",
            "component_mapping": "1000 Genomes Phase 3 EUR genotype-correlation r-squared; trait ABF posterior mass within r2>=0.8 of PBC component leads",
            "classification_boundary": "No-significance results are not interpreted as proof of PBC specificity.",
        },
        "component_leads": COMPONENT_LEADS,
    }
    (PROCESSED_DIR / "10_step3_analysis_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
