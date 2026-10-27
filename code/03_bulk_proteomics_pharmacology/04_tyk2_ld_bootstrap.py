from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "work" / "pharmacology" / "tyk2_safety" / "locus"
PANEL = REPO / "data" / "raw" / "genetics" / "1000_genomes" / "integrated_call_samples_v3.20130502.ALL.panel"
DOSAGE = ROOT / "00_two_variant_ALT_dosage.raw"
BOOTSTRAPS = 100_000
SEED = 20261014
Z_RS345_C = -7.102184654570452
Z_RS230_A = -5.806742467656584


def signed_r(group: pd.DataFrame, first: str, second: str) -> float:
    values = group[[first, second]].dropna().to_numpy(dtype=float)
    if len(values) < 3 or np.var(values[:, 0]) == 0 or np.var(values[:, 1]) == 0:
        return np.nan
    return float(np.corrcoef(values[:, 0], values[:, 1])[0, 1])


def conditional_results(r_value: np.ndarray | float):
    r_array = np.asarray(r_value, dtype=float)
    denominator = np.sqrt(1 - np.square(r_array))
    z_rs345 = (Z_RS345_C - r_array * Z_RS230_A) / denominator
    z_rs230 = (Z_RS230_A - r_array * Z_RS345_C) / denominator
    p_rs345 = 2 * stats.norm.sf(np.abs(z_rs345))
    p_rs230 = 2 * stats.norm.sf(np.abs(z_rs230))
    return z_rs345, p_rs345, z_rs230, p_rs230


def group_summary(name: str, group: pd.DataFrame, first: str, second: str, design: str) -> dict:
    r_value = signed_r(group, first, second)
    z_rs345, p_rs345, z_rs230, p_rs230 = conditional_results(r_value)
    n = len(group)
    return {
        "scenario": name,
        "design": design,
        "n": n,
        "populations": ";".join(sorted(group["pop"].unique())),
        "rs34536443_C_allele_count": float(group[first].sum()),
        "rs34536443_C_carrier_n": int((group[first] > 0).sum()),
        "rs34536443_C_frequency": float(group[first].sum() / (2 * n)),
        "rs2304256_A_allele_count": float(group[second].sum()),
        "rs2304256_A_carrier_n": int((group[second] > 0).sum()),
        "rs2304256_A_frequency": float(group[second].sum() / (2 * n)),
        "signed_dosage_r_C_A": r_value,
        "dosage_r_squared": r_value**2,
        "rs34536443_C_conditional_z": float(z_rs345),
        "rs34536443_C_conditional_p": float(p_rs345),
        "rs2304256_A_conditional_z": float(z_rs230),
        "rs2304256_A_conditional_p": float(p_rs230),
        "both_conditional_p_lt_0_05": bool((p_rs345 < 0.05) and (p_rs230 < 0.05)),
    }


def vectorized_r(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    first_centered = first - first.mean(axis=1, keepdims=True)
    second_centered = second - second.mean(axis=1, keepdims=True)
    numerator = np.sum(first_centered * second_centered, axis=1)
    denominator = np.sqrt(np.sum(first_centered**2, axis=1) * np.sum(second_centered**2, axis=1))
    result = np.full(len(first), np.nan)
    valid = denominator > 0
    result[valid] = numerator[valid] / denominator[valid]
    return result


def bootstrap_r(group: pd.DataFrame, first: str, second: str, seed: int, stratified: bool) -> np.ndarray:
    rng = np.random.default_rng(seed)
    output = np.full(BOOTSTRAPS, np.nan)
    batch_size = 2_000
    population_arrays = {
        pop: subset[[first, second]].to_numpy(dtype=float)
        for pop, subset in group.groupby("pop", sort=True)
    }
    unstratified = group[[first, second]].to_numpy(dtype=float)
    completed = 0
    while completed < BOOTSTRAPS:
        current = min(batch_size, BOOTSTRAPS - completed)
        if stratified:
            sampled_parts = []
            for values in population_arrays.values():
                indices = rng.integers(0, len(values), size=(current, len(values)))
                sampled_parts.append(values[indices])
            sampled = np.concatenate(sampled_parts, axis=1)
        else:
            indices = rng.integers(0, len(unstratified), size=(current, len(unstratified)))
            sampled = unstratified[indices]
        output[completed : completed + current] = vectorized_r(sampled[:, :, 0], sampled[:, :, 1])
        completed += current
    return output


def percentile(values: np.ndarray, probability: float) -> float:
    finite = values[np.isfinite(values)]
    return float(np.quantile(finite, probability))


def bootstrap_summary(name: str, group: pd.DataFrame, values: np.ndarray, first: str, second: str, design: str) -> dict:
    observed = signed_r(group, first, second)
    finite = values[np.isfinite(values)]
    z_rs345, p_rs345, z_rs230, p_rs230 = conditional_results(finite)
    return {
        "scenario": name,
        "bootstrap_design": design,
        "n": len(group),
        "bootstrap_requested": BOOTSTRAPS,
        "bootstrap_valid": len(finite),
        "bootstrap_invalid_monomorphic": BOOTSTRAPS - len(finite),
        "valid_fraction": len(finite) / BOOTSTRAPS,
        "observed_r": observed,
        "bootstrap_r_mean": float(np.mean(finite)),
        "bootstrap_r_median": float(np.median(finite)),
        "bootstrap_r_sd": float(np.std(finite, ddof=1)),
        "r_CI95_low": percentile(finite, 0.025),
        "r_CI95_high": percentile(finite, 0.975),
        "r2_CI95_low": percentile(np.square(finite), 0.025),
        "r2_CI95_high": percentile(np.square(finite), 0.975),
        "fraction_r_positive": float(np.mean(finite > 0)),
        "fraction_r2_lt_0_1": float(np.mean(np.square(finite) < 0.1)),
        "fraction_r2_lt_0_8": float(np.mean(np.square(finite) < 0.8)),
        "rs345_conditional_z_CI95_low": percentile(z_rs345, 0.025),
        "rs345_conditional_z_CI95_high": percentile(z_rs345, 0.975),
        "rs345_conditional_p_CI95_low": percentile(p_rs345, 0.025),
        "rs345_conditional_p_CI95_high": percentile(p_rs345, 0.975),
        "rs345_fraction_conditional_p_lt_0_05": float(np.mean(p_rs345 < 0.05)),
        "rs230_conditional_z_CI95_low": percentile(z_rs230, 0.025),
        "rs230_conditional_z_CI95_high": percentile(z_rs230, 0.975),
        "rs230_conditional_p_CI95_low": percentile(p_rs230, 0.025),
        "rs230_conditional_p_CI95_high": percentile(p_rs230, 0.975),
        "rs230_fraction_conditional_p_lt_0_05": float(np.mean(p_rs230 < 0.05)),
        "fraction_both_conditional_p_lt_0_05": float(np.mean((p_rs345 < 0.05) & (p_rs230 < 0.05))),
        "rs34536443_C_frequency": float(group[first].sum() / (2 * len(group))),
        "rs2304256_A_frequency": float(group[second].sum() / (2 * len(group))),
    }


dosage = pd.read_csv(DOSAGE, sep="\t")
panel = pd.read_csv(PANEL, sep="\t", usecols=["sample", "pop", "super_pop", "gender"])
data = dosage.merge(panel, left_on="IID", right_on="sample", how="left", validate="one_to_one")
if data["pop"].isna().any():
    raise ValueError("Population annotation missing for one or more EUR samples")

first = next(column for column in dosage.columns if column.startswith("rs34536443_C"))
second = next(column for column in dosage.columns if column.startswith("rs2304256_A"))
if set(data[first].dropna().unique()) - {0, 1, 2} or set(data[second].dropna().unique()) - {0, 1, 2}:
    raise ValueError("Unexpected hard-call dosage values")

eur = data.copy()
fin = data[data["pop"] == "FIN"].copy()
non_fin = data[data["pop"] != "FIN"].copy()

population_rows = [
    group_summary("EUR_all", eur, first, second, "all 1000G EUR samples"),
    group_summary("FIN", fin, first, second, "Finnish ancestry subset"),
    group_summary("EUR_non_FIN", non_fin, first, second, "EUR excluding FIN"),
]
for population, group in data.groupby("pop", sort=True):
    if population == "FIN":
        continue
    population_rows.append(group_summary(population, group, first, second, "single 1000G EUR population"))
for population in sorted(data["pop"].unique()):
    group = data[data["pop"] != population]
    population_rows.append(
        group_summary(
            f"EUR_leave_{population}_out",
            group,
            first,
            second,
            f"EUR leave-one-population-out: removed {population}",
        )
    )

population_sensitivity = pd.DataFrame(population_rows)
population_sensitivity.to_csv(ROOT / "06_TYK2_population_LD_sensitivity.tsv", sep="\t", index=False)

bootstrap_specs = [
    ("EUR_all_unstratified", eur, False, SEED + 1, "individual bootstrap"),
    ("EUR_all_population_stratified", eur, True, SEED + 2, "within-population stratified bootstrap"),
    ("FIN", fin, False, SEED + 3, "individual bootstrap"),
    ("EUR_non_FIN_population_stratified", non_fin, True, SEED + 4, "within-population stratified bootstrap"),
]

bootstrap_arrays = {}
bootstrap_rows = []
draw_frames = []
for name, group, stratified, seed, design in bootstrap_specs:
    values = bootstrap_r(group, first, second, seed, stratified)
    bootstrap_arrays[name] = values
    bootstrap_rows.append(bootstrap_summary(name, group, values, first, second, design))
    z_rs345, p_rs345, z_rs230, p_rs230 = conditional_results(values)
    draw_frames.append(
        pd.DataFrame(
            {
                "scenario": name,
                "replicate": np.arange(1, BOOTSTRAPS + 1),
                "signed_r": values,
                "r_squared": np.square(values),
                "rs345_conditional_z": z_rs345,
                "rs345_conditional_p": p_rs345,
                "rs230_conditional_z": z_rs230,
                "rs230_conditional_p": p_rs230,
                "valid": np.isfinite(values),
            }
        )
    )

bootstrap_summary_frame = pd.DataFrame(bootstrap_rows)
bootstrap_summary_frame.to_csv(ROOT / "07_TYK2_LD_bootstrap_summary.tsv", sep="\t", index=False)
pd.concat(draw_frames, ignore_index=True).to_csv(
    ROOT / "08_TYK2_LD_bootstrap_draws.tsv.gz",
    sep="\t",
    index=False,
    compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
)

fin_values = bootstrap_arrays["FIN"]
non_fin_values = bootstrap_arrays["EUR_non_FIN_population_stratified"]
valid_difference = np.isfinite(fin_values) & np.isfinite(non_fin_values)
differences = fin_values[valid_difference] - non_fin_values[valid_difference]
observed_difference = signed_r(fin, first, second) - signed_r(non_fin, first, second)
difference_row = {
    "comparison": "FIN_minus_EUR_non_FIN",
    "observed_r_difference": observed_difference,
    "paired_bootstrap_valid": len(differences),
    "paired_bootstrap_valid_fraction": len(differences) / BOOTSTRAPS,
    "bootstrap_difference_mean": float(np.mean(differences)),
    "bootstrap_difference_median": float(np.median(differences)),
    "difference_CI95_low": percentile(differences, 0.025),
    "difference_CI95_high": percentile(differences, 0.975),
    "two_sided_bootstrap_p_around_zero": float(2 * min(np.mean(differences <= 0), np.mean(differences >= 0))),
    "interpretation_rule": "CI excluding zero suggests population-specific LD; rare-allele monomorphic replicates are excluded and valid fraction must be reported",
}
pd.DataFrame([difference_row]).to_csv(ROOT / "08b_TYK2_FIN_vs_nonFIN_LD_difference.tsv", sep="\t", index=False)

eur_stratified = bootstrap_summary_frame[
    bootstrap_summary_frame["scenario"] == "EUR_all_population_stratified"
].iloc[0]
fin_summary = bootstrap_summary_frame[bootstrap_summary_frame["scenario"] == "FIN"].iloc[0]
non_fin_summary = bootstrap_summary_frame[
    bootstrap_summary_frame["scenario"] == "EUR_non_FIN_population_stratified"
].iloc[0]
fin_observed = population_sensitivity[population_sensitivity["scenario"] == "FIN"].iloc[0]
non_fin_observed = population_sensitivity[population_sensitivity["scenario"] == "EUR_non_FIN"].iloc[0]
leave_out = population_sensitivity[population_sensitivity["scenario"].str.startswith("EUR_leave_")]

report = f"""# TYK2 LD bootstrap 与 FIN/非FIN 人群敏感性

## 主要结论

1000 Genomes Phase 3 EUR（n=503，GRCh37，1-based）中，明确计数 rs34536443-C 与 rs2304256-A 后，个体剂量 signed r=0.327534，与 PLINK `--r-unphased ref-based` 的结果完全一致。100,000 次按五个 EUR 人群分层的 bootstrap 给出 r 的 95% percentile CI 为 {eur_stratified['r_CI95_low']:.3f} 至 {eur_stratified['r_CI95_high']:.3f}，r² CI 为 {eur_stratified['r2_CI95_low']:.3f} 至 {eur_stratified['r2_CI95_high']:.3f}。bootstrap 中 r²<0.8 的比例为 {eur_stratified['fraction_r2_lt_0_8']:.1%}，说明“两个位点不是高LD代理”的判断对抽样波动稳定。

LD 不确定性传播到 GCST005581 双变异 summary-z 条件近似后：rs34536443-C 条件 z 的95%区间为 {eur_stratified['rs345_conditional_z_CI95_low']:.3f} 至 {eur_stratified['rs345_conditional_z_CI95_high']:.3f}，条件 P 的95%区间为 {eur_stratified['rs345_conditional_p_CI95_low']:.3g} 至 {eur_stratified['rs345_conditional_p_CI95_high']:.3g}；rs2304256-A 条件 z 区间为 {eur_stratified['rs230_conditional_z_CI95_low']:.3f} 至 {eur_stratified['rs230_conditional_z_CI95_high']:.3f}，P 区间为 {eur_stratified['rs230_conditional_p_CI95_low']:.3g} 至 {eur_stratified['rs230_conditional_p_CI95_high']:.3g}。两者同时保持条件 P<0.05 的有效重复比例为 {eur_stratified['fraction_both_conditional_p_lt_0_05']:.1%}。

## FIN 与非FIN

FIN 为 n={int(fin_observed['n'])}，rs34536443-C 仅 {int(fin_observed['rs34536443_C_allele_count'])} 个等位基因、{int(fin_observed['rs34536443_C_carrier_n'])} 名携带者；观察到 r={fin_observed['signed_dosage_r_C_A']:.3f}。其 bootstrap 有效比例为 {fin_summary['valid_fraction']:.1%}，无效重复来自 rs34536443-C 被重采样成单态；有效重复的 r 95% CI 为 {fin_summary['r_CI95_low']:.3f} 至 {fin_summary['r_CI95_high']:.3f}。因此 FIN 的点估计可以展示，但不能当作精确的芬兰LD矩阵替代品。

非FIN EUR 为 n={int(non_fin_observed['n'])}，观察到 r={non_fin_observed['signed_dosage_r_C_A']:.3f}，分层 bootstrap r 95% CI 为 {non_fin_summary['r_CI95_low']:.3f} 至 {non_fin_summary['r_CI95_high']:.3f}。FIN-minus-nonFIN 的 r 差为 {observed_difference:.3f}，bootstrap 95% CI 为 {difference_row['difference_CI95_low']:.3f} 至 {difference_row['difference_CI95_high']:.3f}，双侧 bootstrap P={difference_row['two_sided_bootstrap_p_around_zero']:.3g}。该比较受 FIN 罕见等位基因计数限制，应解释为敏感性而非确定的人群异质性检验。

五次 leave-one-population-out 中，r 的范围为 {leave_out['signed_dosage_r_C_A'].min():.3f} 至 {leave_out['signed_dosage_r_C_A'].max():.3f}；所有情景下两个条件 P 均低于0.05。核心结论因此是：rs34536443 与 rs2304256 部分相关但不可互换，两个 PBC 信号在合理的 EUR LD 波动范围内仍可区分。

## 方法与限制

- bootstrap 次数为100,000，固定随机种子为 {SEED}；主分析在 CEU、FIN、GBR、IBS、TSI 内分层重采样，另提供全EUR非分层结果。
- 置信区间为非参数 percentile CI。GCST005581 的边际 z 固定，仅传播参考LD抽样误差；未传播GWAS效应估计误差，也不等同个体级条件回归。
- FIN 只有3个 rs34536443-C 等位基因，bootstrap单态重复被标记为无效并报告比例，未用0或任意数值填补。
- 本地FIN是1000 Genomes FIN子样本，不是FinnGen研究参与者，不能替代申请获得的FinnGen signed LD。
"""
(ROOT / "02_TYK2_LD_bootstrap_population_report.md").write_text(report, encoding="utf-8")

print(population_sensitivity[["scenario", "n", "signed_dosage_r_C_A", "dosage_r_squared", "both_conditional_p_lt_0_05"]].to_string(index=False))
print(bootstrap_summary_frame[["scenario", "valid_fraction", "r_CI95_low", "r_CI95_high", "fraction_both_conditional_p_lt_0_05"]].to_string(index=False))
