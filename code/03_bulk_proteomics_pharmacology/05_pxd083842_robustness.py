from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm


REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "proteomics" / "pxd083842"
ANALYSIS = REPO / "work" / "proteomics" / "pxd083842"
REPORTS = ANALYSIS / "reports"
RNG_SEED = 20261014
PERMUTATIONS = 100_000


def bh_adjust(values: pd.Series) -> pd.Series:
    array = values.to_numpy(dtype=float)
    result = np.full(len(array), np.nan)
    finite = np.isfinite(array)
    ranked_values = array[finite]
    order = np.argsort(ranked_values)
    adjusted = ranked_values[order] * len(ranked_values) / np.arange(1, len(ranked_values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    result[finite] = restored
    return pd.Series(result, index=values.index)


def hedges_g(group_a: np.ndarray, group_b: np.ndarray) -> float:
    n_a, n_b = len(group_a), len(group_b)
    pooled = ((n_a - 1) * np.var(group_a, ddof=1) + (n_b - 1) * np.var(group_b, ddof=1)) / (n_a + n_b - 2)
    if pooled <= 0:
        return 0.0
    correction = 1 - 3 / (4 * (n_a + n_b) - 9)
    return float(correction * (np.mean(group_a) - np.mean(group_b)) / math.sqrt(pooled))


def cliffs_delta(group_a: np.ndarray, group_b: np.ndarray) -> float:
    differences = group_a[:, None] - group_b[None, :]
    return float((np.sum(differences > 0) - np.sum(differences < 0)) / differences.size)


def two_sample(group_a: np.ndarray, group_b: np.ndarray) -> dict:
    welch = stats.ttest_ind(group_a, group_b, equal_var=False)
    mann = stats.mannwhitneyu(group_a, group_b, alternative="two-sided")
    pairwise = group_a[:, None] - group_b[None, :]
    return {
        "n_PBC": len(group_a),
        "n_comparator": len(group_b),
        "mean_PBC": float(np.mean(group_a)),
        "mean_comparator": float(np.mean(group_b)),
        "delta": float(np.mean(group_a) - np.mean(group_b)),
        "median_PBC": float(np.median(group_a)),
        "median_comparator": float(np.median(group_b)),
        "median_difference": float(np.median(group_a) - np.median(group_b)),
        "hodges_lehmann": float(np.median(pairwise)),
        "hedges_g": hedges_g(group_a, group_b),
        "cliffs_delta": cliffs_delta(group_a, group_b),
        "welch_statistic": float(welch.statistic),
        "welch_p": float(welch.pvalue),
        "mann_whitney_U": float(mann.statistic),
        "mann_whitney_p": float(mann.pvalue),
    }


def permutation_p(group_a: np.ndarray, group_b: np.ndarray, seed: int) -> float:
    combined = np.concatenate([group_a, group_b])
    n_a = len(group_a)
    total_sum = float(np.sum(combined))
    observed = abs(float(np.mean(group_a) - np.mean(group_b)))
    rng = np.random.default_rng(seed)
    extreme = 0
    completed = 0
    batch_size = 5_000
    while completed < PERMUTATIONS:
        current = min(batch_size, PERMUTATIONS - completed)
        selected = np.argpartition(rng.random((current, len(combined))), n_a - 1, axis=1)[:, :n_a]
        sum_a = combined[selected].sum(axis=1)
        differences = sum_a / n_a - (total_sum - sum_a) / (len(combined) - n_a)
        extreme += int(np.sum(np.abs(differences) >= observed - 1e-12))
        completed += current
    return float((extreme + 1) / (PERMUTATIONS + 1))


def regression_results(group_a: np.ndarray, group_b: np.ndarray) -> list[dict]:
    response = np.concatenate([group_b, group_a])
    indicator = np.concatenate([np.zeros(len(group_b)), np.ones(len(group_a))])
    design = sm.add_constant(indicator)
    ols = sm.OLS(response, design).fit(cov_type="HC3")
    rlm = sm.RLM(response, design, M=sm.robust.norms.HuberT()).fit()
    return [
        {
            "method": "OLS_HC3",
            "estimate": float(ols.params[1]),
            "standard_error": float(ols.bse[1]),
            "statistic": float(ols.tvalues[1]),
            "p_value": float(ols.pvalues[1]),
        },
        {
            "method": "Huber_RLM",
            "estimate": float(rlm.params[1]),
            "standard_error": float(rlm.bse[1]),
            "statistic": float(rlm.tvalues[1]),
            "p_value": float(rlm.pvalues[1]),
        },
    ]


def cap_values(values: np.ndarray, method: str) -> np.ndarray:
    if method == "winsor_5_95":
        lower, upper = np.quantile(values, [0.05, 0.95])
    else:
        median = np.median(values)
        mad = np.median(np.abs(values - median))
        scale = 1.4826 * mad
        lower, upper = median - 3 * scale, median + 3 * scale
    return np.clip(values, lower, upper)


def read_pxd_matrix(path: Path):
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str)
    sample_ids = raw.iloc[0, 2:].tolist()
    groups = pd.Series(raw.iloc[1, 2:].tolist(), index=sample_ids, name="group")
    subjects = pd.Series(raw.iloc[2, 2:].tolist(), index=sample_ids, name="subject_id")
    genes = raw.iloc[3:, 1].astype(str).str.split("~", n=1).str[0]
    values = raw.iloc[3:, 2:].apply(pd.to_numeric, errors="coerce")
    values.columns = sample_ids
    values.index = genes
    return values.groupby(level=0).median(), groups, subjects


def reactome_genes(path: Path) -> list[str]:
    participants = json.loads(path.read_text(encoding="utf-8-sig"))
    genes = set()
    for participant in participants:
        for entity in participant.get("refEntities", []):
            if entity.get("schemaClass") not in {"ReferenceGeneProduct", "ReferenceIsoform"}:
                continue
            tokens = str(entity.get("displayName", "")).split()
            if len(tokens) >= 2:
                genes.add(tokens[-1])
    return sorted(genes)


def comparator_samples(groups: pd.Series, comparator: str) -> pd.Index:
    if comparator == "NonPBC":
        return groups.index[groups != "PBC"]
    return groups.index[groups == comparator]


def make_pathway_scores(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    selected = matrix.loc[genes]
    standardized = selected.sub(selected.mean(axis=1), axis=0).div(selected.std(axis=1).replace(0, np.nan), axis=0)
    return standardized.mean(axis=0, skipna=True)


liver_path = RAW / "PXD083842" / "liver_protein_data_minprob_processed.tsv.txt"
plasma_path = RAW / "PXD083842" / "plasma_protein_data_minprob_processed.tsv.txt"
liver, liver_groups, liver_subjects = read_pxd_matrix(liver_path)
plasma, plasma_groups, plasma_subjects = read_pxd_matrix(plasma_path)

comparators = ["Healthy", "NAFLD", "NonPBC"]
tyk2 = liver.loc["TYK2"].astype(float)

leave_one_out_rows = []
for comparator in comparators:
    selected_samples = liver_groups.index[(liver_groups == "PBC") | liver_groups.index.isin(comparator_samples(liver_groups, comparator))]
    for dropped_sample in selected_samples:
        retained = selected_samples[selected_samples != dropped_sample]
        group_a = tyk2[retained[liver_groups[retained] == "PBC"]].to_numpy()
        group_b = tyk2[retained[liver_groups[retained] != "PBC"]].to_numpy()
        row = {
            "dataset": "PXD083842",
            "tissue": "liver",
            "feature": "TYK2",
            "contrast": f"PBC_vs_{comparator}",
            "dropped_sample": dropped_sample,
            "dropped_group": liver_groups[dropped_sample],
            "dropped_subject": liver_subjects[dropped_sample],
        }
        row.update(two_sample(group_a, group_b))
        leave_one_out_rows.append(row)

leave_one_out = pd.DataFrame(leave_one_out_rows)
leave_one_out.to_csv(ANALYSIS / "10_PXD083842_TYK2_leave_one_out.tsv", sep="\t", index=False)

robust_rows = []
for contrast_index, comparator in enumerate(comparators):
    group_a = tyk2[liver_groups.index[liver_groups == "PBC"]].to_numpy()
    group_b = tyk2[comparator_samples(liver_groups, comparator)].to_numpy()
    basic = two_sample(group_a, group_b)
    common = {
        "dataset": "PXD083842",
        "tissue": "liver",
        "feature": "TYK2",
        "contrast": f"PBC_vs_{comparator}",
        "n_PBC": len(group_a),
        "n_comparator": len(group_b),
    }
    for method, estimate, statistic, p_value in [
        ("Welch", basic["delta"], basic["welch_statistic"], basic["welch_p"]),
        ("Mann_Whitney", basic["hodges_lehmann"], basic["mann_whitney_U"], basic["mann_whitney_p"]),
        ("Permutation_100000", basic["delta"], np.nan, permutation_p(group_a, group_b, RNG_SEED + contrast_index)),
        ("Hodges_Lehmann", basic["hodges_lehmann"], basic["cliffs_delta"], basic["mann_whitney_p"]),
    ]:
        robust_rows.append({**common, "method": method, "estimate": estimate, "standard_error": np.nan, "statistic": statistic, "p_value": p_value})
    for regression in regression_results(group_a, group_b):
        robust_rows.append({**common, **regression})
    pooled = np.concatenate([group_a, group_b])
    for cap_index, cap_method in enumerate(["winsor_5_95", "MAD_3sigma"]):
        capped = cap_values(pooled, cap_method)
        capped_a = capped[: len(group_a)]
        capped_b = capped[len(group_a) :]
        capped_stats = two_sample(capped_a, capped_b)
        robust_rows.append(
            {
                **common,
                "method": cap_method,
                "estimate": capped_stats["delta"],
                "standard_error": np.nan,
                "statistic": capped_stats["welch_statistic"],
                "p_value": capped_stats["welch_p"],
                "secondary_nonparametric_p": capped_stats["mann_whitney_p"],
            }
        )

robust_methods = pd.DataFrame(robust_rows)
robust_methods["FDR_within_method_family"] = robust_methods.groupby("method")["p_value"].transform(bh_adjust)
robust_methods.to_csv(ANALYSIS / "11_PXD083842_TYK2_robust_methods.tsv", sep="\t", index=False)

ifng_reference = reactome_genes(RAW / "reactome_participants_R-HSA-877300.json")
ifng_genes = sorted(gene for gene in ifng_reference if gene in liver.index and not gene.startswith("HLA-") and gene != "B2M")

leave_gene_rows = []
for omitted_gene in ifng_genes:
    retained_genes = [gene for gene in ifng_genes if gene != omitted_gene]
    scores = make_pathway_scores(liver, retained_genes)
    for comparator in comparators:
        group_a = scores[liver_groups.index[liver_groups == "PBC"]].to_numpy()
        group_b = scores[comparator_samples(liver_groups, comparator)].to_numpy()
        row = {
            "dataset": "PXD083842",
            "tissue": "liver",
            "pathway": "Reactome_IFNG_without_HLA_B2M",
            "contrast": f"PBC_vs_{comparator}",
            "omitted_gene": omitted_gene,
            "genes_retained": len(retained_genes),
        }
        row.update(two_sample(group_a, group_b))
        leave_gene_rows.append(row)

leave_gene = pd.DataFrame(leave_gene_rows)
leave_gene["welch_FDR_across_omitted_genes"] = leave_gene.groupby("contrast")["welch_p"].transform(bh_adjust)
leave_gene.to_csv(ANALYSIS / "12_PXD083842_IFNG_leave_one_gene_out.tsv", sep="\t", index=False)

leave_sample_rows = []
for comparator in comparators:
    selected_samples = liver_groups.index[(liver_groups == "PBC") | liver_groups.index.isin(comparator_samples(liver_groups, comparator))]
    for dropped_sample in selected_samples:
        retained_all = liver.columns[liver.columns != dropped_sample]
        recalculated_scores = make_pathway_scores(liver[retained_all], ifng_genes)
        group_a = recalculated_scores[liver_groups.index[(liver_groups == "PBC") & (liver_groups.index != dropped_sample)]].to_numpy()
        group_b = recalculated_scores[comparator_samples(liver_groups, comparator).difference([dropped_sample], sort=False)].to_numpy()
        row = {
            "dataset": "PXD083842",
            "tissue": "liver",
            "pathway": "Reactome_IFNG_without_HLA_B2M",
            "contrast": f"PBC_vs_{comparator}",
            "dropped_sample": dropped_sample,
            "dropped_group": liver_groups[dropped_sample],
            "dropped_subject": liver_subjects[dropped_sample],
            "genes_measured": len(ifng_genes),
        }
        row.update(two_sample(group_a, group_b))
        leave_sample_rows.append(row)

leave_sample = pd.DataFrame(leave_sample_rows)
leave_sample.to_csv(ANALYSIS / "13_PXD083842_IFNG_leave_one_sample_out.tsv", sep="\t", index=False)

subject_rows = []
all_subjects = sorted(set(liver_subjects) | set(plasma_subjects))
for subject in all_subjects:
    liver_samples = liver_subjects.index[liver_subjects == subject].tolist()
    plasma_samples = plasma_subjects.index[plasma_subjects == subject].tolist()
    liver_sample = liver_samples[0] if liver_samples else None
    plasma_sample = plasma_samples[0] if plasma_samples else None
    subject_rows.append(
        {
            "subject_id": subject,
            "group": liver_groups[liver_sample] if liver_sample else plasma_groups[plasma_sample],
            "liver_sample": liver_sample,
            "plasma_sample": plasma_sample,
            "both_tissues_available": bool(liver_sample and plasma_sample),
            "liver_TYK2_measured": bool(liver_sample and "TYK2" in liver.index),
            "plasma_TYK2_measured": bool(plasma_sample and "TYK2" in plasma.index),
            "liver_TYK2": float(liver.loc["TYK2", liver_sample]) if liver_sample and "TYK2" in liver.index else np.nan,
            "plasma_TYK2": float(plasma.loc["TYK2", plasma_sample]) if plasma_sample and "TYK2" in plasma.index else np.nan,
        }
    )

subject_overlap = pd.DataFrame(subject_rows)
subject_overlap.to_csv(ANALYSIS / "14_PXD083842_liver_plasma_subject_overlap.tsv", sep="\t", index=False)

summary_rows = []
for contrast in leave_one_out["contrast"].unique():
    subset = leave_one_out[leave_one_out["contrast"] == contrast]
    summary_rows.append(
        {
            "analysis": "TYK2_leave_one_sample_out",
            "contrast": contrast,
            "iterations": len(subset),
            "min_delta": subset["delta"].min(),
            "max_delta": subset["delta"].max(),
            "positive_direction_fraction": (subset["delta"] > 0).mean(),
            "min_hedges_g": subset["hedges_g"].min(),
            "max_welch_p": subset["welch_p"].max(),
            "fraction_welch_p_lt_0_05": (subset["welch_p"] < 0.05).mean(),
            "max_mann_whitney_p": subset["mann_whitney_p"].max(),
            "fraction_mann_whitney_p_lt_0_05": (subset["mann_whitney_p"] < 0.05).mean(),
        }
    )
for contrast in leave_gene["contrast"].unique():
    subset = leave_gene[leave_gene["contrast"] == contrast]
    summary_rows.append(
        {
            "analysis": "IFNG_leave_one_gene_out",
            "contrast": contrast,
            "iterations": len(subset),
            "min_delta": subset["delta"].min(),
            "max_delta": subset["delta"].max(),
            "positive_direction_fraction": (subset["delta"] > 0).mean(),
            "min_hedges_g": subset["hedges_g"].min(),
            "max_welch_p": subset["welch_p"].max(),
            "fraction_welch_p_lt_0_05": (subset["welch_p"] < 0.05).mean(),
            "max_mann_whitney_p": subset["mann_whitney_p"].max(),
            "fraction_mann_whitney_p_lt_0_05": (subset["mann_whitney_p"] < 0.05).mean(),
        }
    )
for contrast in leave_sample["contrast"].unique():
    subset = leave_sample[leave_sample["contrast"] == contrast]
    summary_rows.append(
        {
            "analysis": "IFNG_leave_one_sample_out",
            "contrast": contrast,
            "iterations": len(subset),
            "min_delta": subset["delta"].min(),
            "max_delta": subset["delta"].max(),
            "positive_direction_fraction": (subset["delta"] > 0).mean(),
            "min_hedges_g": subset["hedges_g"].min(),
            "max_welch_p": subset["welch_p"].max(),
            "fraction_welch_p_lt_0_05": (subset["welch_p"] < 0.05).mean(),
            "max_mann_whitney_p": subset["mann_whitney_p"].max(),
            "fraction_mann_whitney_p_lt_0_05": (subset["mann_whitney_p"] < 0.05).mean(),
        }
    )

robustness_summary = pd.DataFrame(summary_rows)
robustness_summary.to_csv(ANALYSIS / "15_PXD083842_robustness_summary.tsv", sep="\t", index=False)

tyk2_healthy = robust_methods[robust_methods["contrast"] == "PBC_vs_Healthy"].set_index("method")
tyk2_nafld = robust_methods[robust_methods["contrast"] == "PBC_vs_NAFLD"].set_index("method")
tyk2_loo_healthy = robustness_summary[(robustness_summary["analysis"] == "TYK2_leave_one_sample_out") & (robustness_summary["contrast"] == "PBC_vs_Healthy")].iloc[0]
ifng_gene_healthy = robustness_summary[(robustness_summary["analysis"] == "IFNG_leave_one_gene_out") & (robustness_summary["contrast"] == "PBC_vs_Healthy")].iloc[0]
ifng_sample_healthy = robustness_summary[(robustness_summary["analysis"] == "IFNG_leave_one_sample_out") & (robustness_summary["contrast"] == "PBC_vs_Healthy")].iloc[0]
ifng_gene_worst_healthy = leave_gene[leave_gene["contrast"] == "PBC_vs_Healthy"].sort_values("welch_p", ascending=False).iloc[0]
overlap_count = int(subject_overlap["both_tissues_available"].sum())

report = f"""# PXD083842 稳健性分析

## 结论

PXD083842 的肝组织 TYK2 升高不是由单个样本或单一种参数检验制造。PBC 对 healthy 的原始均值差为 {tyk2_healthy.loc['Welch', 'estimate']:.3f} log2；100,000 次标签置换 P={tyk2_healthy.loc['Permutation_100000', 'p_value']:.4g}，HC3 稳健标准误 P={tyk2_healthy.loc['OLS_HC3', 'p_value']:.4g}，Huber 稳健回归 P={tyk2_healthy.loc['Huber_RLM', 'p_value']:.4g}，Mann–Whitney P={tyk2_healthy.loc['Mann_Whitney', 'p_value']:.4g}。5/95% winsorization 后 P={tyk2_healthy.loc['winsor_5_95', 'p_value']:.4g}。

逐样本删除后，TYK2 的 PBC-minus-healthy 差异范围为 {tyk2_loo_healthy['min_delta']:.3f} 至 {tyk2_loo_healthy['max_delta']:.3f}，{tyk2_loo_healthy['positive_direction_fraction']:.1%} 保持正方向；Welch P<0.05 的比例为 {tyk2_loo_healthy['fraction_welch_p_lt_0_05']:.1%}，最差 P={tyk2_loo_healthy['max_welch_p']:.4g}。因此 TYK2 可保留为主文的外部肝蛋白证据，但它仍是疾病状态关联，不是药物靶点因果效应。

TYK2 的 PBC 对 MASLD 差异仅为 {tyk2_nafld.loc['Welch', 'estimate']:.3f} log2，Welch P={tyk2_nafld.loc['Welch', 'p_value']:.3g}，置换 P={tyk2_nafld.loc['Permutation_100000', 'p_value']:.3g}，各稳健方法均未达到 0.05。故该结果支持 PBC 相对健康肝的升高，不能声称 TYK2 是区分 PBC 与其他慢性脂肪性肝病的特异蛋白标志物。

去除 HLA/B2M 的 IFN-γ 程序包含 {len(ifng_genes)} 个实测蛋白。逐基因删除后，PBC 对 healthy 的方向保持率为 {ifng_gene_healthy['positive_direction_fraction']:.1%}，Welch P<0.05 比例为 {ifng_gene_healthy['fraction_welch_p_lt_0_05']:.1%}；删除最具影响的 {ifng_gene_worst_healthy['omitted_gene']} 后仍为 P={ifng_gene_worst_healthy['welch_p']:.4g}。逐样本删除后的方向保持率为 {ifng_sample_healthy['positive_direction_fraction']:.1%}，Welch P<0.05 比例为 {ifng_sample_healthy['fraction_welch_p_lt_0_05']:.1%}，最差 P={ifng_sample_healthy['max_welch_p']:.4g}。这检验了结果是否由单一 HLA/B2M、单个通路成员或单个样本主导。

## 配对与数据限制

- 肝和血浆样本按公开 subject ID 匹配后共有 {overlap_count} 位同时具备两种组织；公开矩阵不是全队列严格配对。
- 血浆矩阵中 TYK2 是否可定量已在 `14_PXD083842_liver_plasma_subject_overlap.tsv` 逐人记录；不能把肝 TYK2 结果描述为跨组织复制。
- 两个矩阵是提交者发布的 MinProb-imputed 处理矩阵。公开数据不含可重建完整缺失模式的未填补定量矩阵，因此本轮无法诚实完成“无填补复算”。
- 项目状态仍为 PARTIAL，且截至本次审计没有关联同行评议论文。建议主文中把它写作独立公开蛋白组的稳健外部支持，并在局限性中明确新发布、处理后矩阵和 MinProb 填补。

## 统计口径

- 主要效应同时报告均值差、Hedges g、Hodges–Lehmann 差和 Cliff's delta。
- 参数检验使用 Welch t；非参数检验使用双侧 Mann–Whitney；稳健性使用 HC3 OLS、Huber RLM、5/95% winsorization、3-MAD 限幅和 100,000 次固定随机种子标签置换。
- 原始全蛋白差异的 FDR 仍以既有全蛋白 BH 校正为主；稳健方法是同一预设 TYK2 假设的敏感性分析，不能被当作多次独立复制。
"""
(REPORTS / "02_PXD083842_robustness_report.md").write_text(report, encoding="utf-8")

print(f"TYK2 leave-one-out rows: {len(leave_one_out)}")
print(f"IFNG genes: {len(ifng_genes)}; gene-LOO rows: {len(leave_gene)}; sample-LOO rows: {len(leave_sample)}")
print(f"Matched liver/plasma subjects: {overlap_count}")
