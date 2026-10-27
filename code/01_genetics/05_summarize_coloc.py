from __future__ import annotations

from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[2]
RESULT = REPO / "work" / "genetics" / "discovery" / "results"


def best_signal_pair(
    frame: pd.DataFrame, group_columns: list[str] | None = None
) -> pd.DataFrame:
    if group_columns is None:
        group_columns = ["analysis_id"]
    return (
        frame.sort_values(group_columns + ["PP.H4"], ascending=[True] * len(group_columns) + [False])
        .groupby(group_columns, as_index=False)
        .first()
    )


def main() -> None:
    pbc = pd.read_csv(RESULT / "19_pbc_susie_multisignal_finemapping.tsv", sep="\t")
    pbc_cs = pd.read_csv(RESULT / "19b_pbc_susie_credible_set_details.tsv", sep="\t")
    pair_qc = pd.read_csv(RESULT / "20_onek1k_susie_concordance.tsv", sep="\t")
    coloc = pd.read_csv(RESULT / "21_coloc_susie_signal_pairs.tsv", sep="\t")
    sensitivity = pd.read_csv(RESULT / "22_coloc_susie_prior_sensitivity.tsv", sep="\t")
    review_qc = pd.read_csv(
        RESULT / "20_onek1k_susie_concordance_l3_review.tsv", sep="\t"
    )
    review_coloc = pd.read_csv(
        RESULT / "21_coloc_susie_signal_pairs_l3_review.tsv", sep="\t"
    )
    review_sensitivity = pd.read_csv(
        RESULT / "22_coloc_susie_prior_sensitivity_l3_review.tsv", sep="\t"
    )

    strict = best_signal_pair(coloc)
    review_strict = best_signal_pair(review_coloc)
    sa09_review = review_strict.loc[review_strict["analysis_id"].eq("SA09")].copy()
    strict = strict.loc[~strict["analysis_id"].eq("SA09")]
    strict = pd.concat([strict, sa09_review], ignore_index=True)

    primary_qc = pair_qc.copy()
    review_qc_indexed = review_qc.set_index("analysis_id")
    for analysis_id in ["SA09", "SA19"]:
        if analysis_id in review_qc_indexed.index:
            replacement = review_qc_indexed.loc[analysis_id]
            mask = primary_qc["analysis_id"].eq(analysis_id)
            for column in [
                "pbc_converged",
                "pbc_credible_set_n",
                "eqtl_converged",
                "eqtl_credible_set_n",
            ]:
                primary_qc.loc[mask, column] = replacement[column]

    best_by_prior = best_signal_pair(sensitivity, ["analysis_id", "p12"])
    review_best_by_prior = best_signal_pair(review_sensitivity, ["analysis_id", "p12"])
    best_by_prior = best_by_prior.loc[~best_by_prior["analysis_id"].eq("SA09")]
    best_by_prior = pd.concat(
        [best_by_prior, review_best_by_prior.loc[review_best_by_prior["analysis_id"].eq("SA09")]],
        ignore_index=True,
    )
    prior_wide = best_by_prior.pivot(index="analysis_id", columns="p12", values="PP.H4")
    prior_wide = prior_wide.rename(
        columns={1e-6: "best_pph4_p12_1e6", 5e-6: "best_pph4_p12_5e6", 1e-5: "best_pph4_p12_1e5"}
    ).reset_index()

    summary = strict.merge(
        primary_qc[
            [
                "analysis_id",
                "pbc_converged",
                "eqtl_converged",
                "pbc_credible_set_n",
                "eqtl_credible_set_n",
                "pbc_rss_s_hat",
                "eqtl_rss_s_hat",
                "pbc_kriging_outlier_n",
                "eqtl_kriging_outlier_n",
            ]
        ],
        on="analysis_id",
        how="left",
    ).merge(prior_wide, on="analysis_id", how="left")
    summary["model_source"] = "L10_main"
    summary.loc[summary["analysis_id"].eq("SA09"), "model_source"] = "L3_converged_review"
    summary["final_evidence_class"] = summary["evidence_class"]
    summary.loc[summary["analysis_id"].eq("SA19"), "final_evidence_class"] = (
        "excluded_nonconvergent_eqtl"
    )
    summary["strict_prior_robust"] = summary["evidence_class"].isin(["strong", "candidate"])
    summary["interpretation_note"] = "usable"
    summary.loc[summary["analysis_id"].eq("SA09"), "interpretation_note"] = (
        "L10 eQTL did not converge; L3 converged review used"
    )
    summary.loc[summary["analysis_id"].eq("SA19"), "interpretation_note"] = (
        "eQTL SuSiE failed to converge at L10 and L3; excluded from main inference"
    )
    summary = summary.sort_values(
        ["final_evidence_class", "PP.H4"], ascending=[True, False]
    )
    output_columns = [
        "analysis_id",
        "pbc_locus",
        "gene_symbol",
        "gene_id",
        "sample_group",
        "sample_size",
        "variant_n",
        "model_source",
        "pbc_signal",
        "eqtl_signal",
        "pbc_converged",
        "eqtl_converged",
        "pbc_credible_set_n",
        "eqtl_credible_set_n",
        "PP.H3",
        "PP.H4",
        "shared_fraction_H4_over_H3_H4",
        "top_pbc_snp",
        "top_eqtl_snp",
        "top_shared_snp",
        "top_shared_snp_posterior",
        "pbc_credible_set_size",
        "eqtl_credible_set_size",
        "credible_set_overlap_n",
        "credible_set_jaccard",
        "best_pph4_p12_1e6",
        "best_pph4_p12_5e6",
        "best_pph4_p12_1e5",
        "evidence_class",
        "final_evidence_class",
        "interpretation_note",
    ]
    summary[output_columns].to_csv(
        RESULT / "23a_coloc_susie_combination_summary.tsv", sep="\t", index=False
    )

    pbc_l10 = pbc.loc[pbc["model_L"].eq(10)].copy()
    cs_l10 = pbc_cs.loc[pbc_cs["model_L"].eq(10)].copy()
    pbc_l10["credible_set_sizes"] = pbc_l10["locus_id"].map(
        cs_l10.groupby("locus_id")["credible_set_size"]
        .apply(lambda values: ";".join(map(str, values.astype(int))))
    )
    pbc_l10["minimum_cs_purity"] = pbc_l10["locus_id"].map(
        cs_l10.groupby("locus_id")["minimum_abs_correlation"].min()
    )
    pbc_l10.to_csv(RESULT / "23b_pbc_locus_finemapping_summary.tsv", sep="\t", index=False)

    class_counts = summary["final_evidence_class"].value_counts().to_dict()
    pbc_rows = []
    for row in pbc_l10.sort_values("locus_id").itertuples():
        pbc_rows.append(
            f"| {row.pbc_locus} | {int(row.credible_set_n)} | {row.credible_set_sizes} | "
            f"{row.top_snp} | {row.top_pip:.3f} | {row.rss_s_hat:.3f} | "
            f"{int(row.kriging_outlier_n)} |"
        )

    final_rows = []
    rank_order = {
        "strong": 0,
        "candidate": 1,
        "exploratory": 2,
        "no_support": 3,
        "excluded_nonconvergent_eqtl": 4,
    }
    ranked = summary.assign(
        rank=summary["final_evidence_class"].map(rank_order).fillna(9)
    ).sort_values(["rank", "PP.H4"], ascending=[True, False])
    for _, row in ranked.iterrows():
        final_rows.append(
            f"| {row['analysis_id']} | {row['pbc_locus']} | {row['gene_symbol']} | {row['sample_group']} | "
            f"{row['PP.H4']:.3f} | {row['shared_fraction_H4_over_H3_H4']:.3f} | "
            f"{int(row['credible_set_overlap_n'])} | {row['final_evidence_class']} |"
        )

    report = f"""# PBC 局部 LD、多信号精细定位与 coloc.susie 验证报告

## 决策结论

**继续 PBC 主线（GO），不换题。** 修正 LD 计数等位基因方向后，10 个 PBC 位点的 SuSiE 模型全部收敛，LD–汇总统计一致性诊断恢复到可接受范围；23 个基因–细胞组合中，22 个可形成可解释结果，NFKB1/CD4-TCM 因 eQTL SuSiE 在 `L=10` 与 `L=3` 均不收敛而排除。

严格先验 `p12=1e-6` 下，最终得到：

- 强证据组合 {class_counts.get('strong', 0)} 个；
- 候选组合 {class_counts.get('candidate', 0)} 个；
- 探索性组合 {class_counts.get('exploratory', 0)} 个；
- 无支持组合 {class_counts.get('no_support', 0)} 个；
- 模型不稳定并排除 {class_counts.get('excluded_nonconvergent_eqtl', 0)} 个。

这些是“组合数”，不是独立遗传证据数。同一位点的不同基因或细胞类型必须合并为一个遗传单元解释。

## 关键质量控制

- PBC GWAS 与 1000G 均为 GRCh37、1-based；OneK1K 为 GRCh38、1-based。跨版本只以 rsID 加等位基因连接，不直接比较坐标。
- 23 个组合严格对齐 1,147–3,763 个 SNP，共 58,042 条组合–变异记录；10 个 PBC 位点级输入包含 1,586–4,054 个 SNP。
- 所有效应统一到 1000G PVAR ALT；LD 使用 PLINK2 `--r-unphased ... ref-based`。REF–REF 与 ALT–ALT 的相关符号相同。
- 初始试跑使用 PLINK2 默认“主等位基因”相关，造成不同 SNP 间 LD 符号切换，表现为 `s_hat=0.41–0.90`、大量 PIP=1 和普遍不收敛。该错误已定位并完全重算，错误试跑结果未用于结论。
- 修正后 20/20 个 PBC `L=3/L=10` 模型收敛；位点级 `s_hat` 为 {pbc_l10.rss_s_hat.min():.3f}–{pbc_l10.rss_s_hat.max():.3f}，kriging 异常变异为 0–1 个。
- LD 矩阵 23/23 通过尺寸、SNP 集合、顺序、有限值、对角线和对称性检查。float32 读取后显式对称化；未通过大 ridge 掩盖近奇异结构。

## PBC 多信号精细定位

| PBC 位点 | L=10 可信集数 | 可信集大小 | 最高 PIP SNP | 最高 PIP | s_hat | kriging 异常数 |
|---|---:|---|---|---:|---:|---:|
{chr(10).join(pbc_rows)}

主要发现：

- `1p31.3/IL12RB2` 稳定检出 3 个 PBC 信号，`L=3` 与 `L=10` 一致；三个 95% 可信集分别为 2、3、2 个 SNP，纯度高。
- 其余 9 个位点在 `L=3` 与 `L=10` 下均为单信号，说明此前单信号 coloc 在这些位点的主要局限不是 PBC 侧多重独立信号。
- `5q31.3/NDFIP1`、`4q24` 和 `7q34` 的 PBC 可信集较宽，最高 PIP 较低，定位精度仍受外部 LD 和样本量限制。

## coloc.susie 组合级结果

| ID | 位点 | 基因 | 细胞 | PP.H4 | H4/(H3+H4) | CS 重叠 SNP | 最终分级 |
|---|---|---|---|---:|---:|---:|---|
{chr(10).join(final_rows)}

## 审稿人视角下的主结论

1. **IL12RB2/NK**：PBC 位点为真实多信号结构。`L=10` eQTL 模型未收敛，但 `L=3` 收敛复核支持 PBC signal 1 与 eQTL signal 1 共享，严格先验 PP.H4=0.971；可列为强证据，但方法中必须预先声明采用收敛的 `L=3` 复核结果。
2. **LINC00299/ID2 竞争**：LINC00299/NK 为强证据（PP.H4=0.959），ID2/NK 无支持（PP.H4<0.001）。在当前 OneK1K 数据中，遗传证据明确偏向 LINC00299，而不是 ID2。
3. **FCRL3/FCRL3-AS1**：四个组合均为强证据，但共享同一 PBC 遗传信号，不能计为四条独立证据，也不能仅凭共定位区分 FCRL3 与反义转录本。
4. **ELMO1**：两个 CD4 亚群均为严格先验候选，较宽松先验转为强证据；属于先验敏感但跨细胞复现的候选。
5. **ORMDL3/GSDMB**：两个基因均为候选且共享同一位点，现阶段不能声称 ORMDL3 优于 GSDMB；IKZF3 因 OneK1K 无可用关联数据，仍不可比较。
6. **CD226/DOK6**：四个组合均为候选且可信集重叠明显，但同位点基因竞争未解决，应以“CD226–DOK6 调控单元”表述。
7. **NDFIP1**：两个细胞类型均为候选，可信集 Jaccard=0.90，跨细胞一致，但 PBC 可信集较宽。
8. **4q24**：NFKB1 eQTL SuSiE 在 `L=10` 与 `L=3` 均未收敛，排除；MANBA 仅探索性。该位点暂不进入核心排名。
9. **11q13.1**：AP003774.1 为候选/探索性，CCDC88B 为探索性，证据倾向 AP003774.1，但不足以完成基因归属。

## 主线排序建议

按“独立遗传单元”而不是组合数排序：

1. 第一梯队：`IL12RB2`、`LINC00299`、`FCRL3/FCRL3-AS1`。
2. 第二梯队：`ELMO1`、`ORMDL3/GSDMB`、`CD226/DOK6`、`NDFIP1`、`ZC3HAV1`。
3. 第三梯队：`AP003774.1/CCDC88B`、`MANBA`。
4. 暂时排除：`ID2`（当前数据不支持）、`NFKB1`（模型不收敛）、`IKZF3`（OneK1K 无可用数据）。

## 方法边界与下一步

- `coloc.susie` 支持共享遗传信号，不等于已证实因果基因或药物靶点。
- 1000G EUR 仅 503 人，仍是外部 LD；核心位点后续应争取 PBC 同队列 LD，或至少用 UKB/更大欧洲参考面板复核。
- 下一阶段不应继续堆叠同类共定位，而应对第一、二梯队开展：独立 eQTL 数据复现、组织/细胞特异性验证、方向一致性、蛋白/药物可成药性与安全性整合。
- 论文主叙事可以保留“单细胞 eQTL 指导的遗传支持靶点优先排序”，但结果单位必须是位点级遗传单元，并将同位点竞争基因作为不确定性显式展示。
"""
    (RESULT / "23_pbc_multisignal_validation_report.md").write_text(report, encoding="utf-8")
    print(summary["final_evidence_class"].value_counts().to_string())
    print(RESULT / "23_pbc_multisignal_validation_report.md")


if __name__ == "__main__":
    main()
