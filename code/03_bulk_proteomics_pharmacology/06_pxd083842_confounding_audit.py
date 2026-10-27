from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pypdf import PdfReader
from scipy import stats


REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw" / "proteomics" / "pxd083842"
SUPPLEMENTS = RAW / "nature_supplements"
ANALYSIS = REPO / "work" / "proteomics" / "pxd083842"
REPORTS = ANALYSIS / "reports"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_matrix_metadata(path: Path) -> tuple[pd.Series, pd.Series]:
    frame = pd.read_csv(path, sep="\t", header=None, dtype=str)
    sample_ids = frame.iloc[0, 2:].tolist()
    groups = pd.Series(frame.iloc[1, 2:].tolist(), index=sample_ids, name="group")
    subjects = pd.Series(frame.iloc[2, 2:].tolist(), index=sample_ids, name="subject_id")
    return groups, subjects


def read_matrix(path: Path) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    frame = pd.read_csv(path, sep="\t", header=None, dtype=str)
    sample_ids = frame.iloc[0, 2:].tolist()
    groups = pd.Series(frame.iloc[1, 2:].tolist(), index=sample_ids, name="group")
    subjects = pd.Series(frame.iloc[2, 2:].tolist(), index=sample_ids, name="subject_id")
    genes = frame.iloc[3:, 1].astype(str).str.split("~", n=1).str[0]
    values = frame.iloc[3:, 2:].apply(pd.to_numeric, errors="coerce")
    values.columns = sample_ids
    values.index = genes
    return values.groupby(level=0).median(), groups, subjects


def significant_sheet(path: Path, sheet: str) -> pd.DataFrame:
    frame = pd.read_excel(path, sheet_name=sheet, header=2)
    frame = frame.dropna(how="all")
    if "GeneName" in frame.columns:
        frame["GeneName"] = frame["GeneName"].astype(str)
    return frame


def first_gene_row(path: Path, sheets: list[str], gene: str) -> tuple[pd.Series | None, str | None]:
    for sheet in sheets:
        frame = significant_sheet(path, sheet)
        if "GeneName" not in frame.columns:
            continue
        selected = frame[frame["GeneName"].str.upper() == gene.upper()]
        if not selected.empty:
            return selected.iloc[0], sheet
    return None, None


def sheet_gene_sets(path: Path, up_sheet: str, down_sheet: str) -> tuple[set[str], set[str]]:
    up = significant_sheet(path, up_sheet)
    down = significant_sheet(path, down_sheet)
    up_genes = set(up.get("GeneName", pd.Series(dtype=str)).dropna().astype(str))
    down_genes = set(down.get("GeneName", pd.Series(dtype=str)).dropna().astype(str))
    return up_genes, down_genes


def reactome_genes(path: Path) -> set[str]:
    participants = json.loads(path.read_text(encoding="utf-8-sig"))
    genes: set[str] = set()
    for participant in participants:
        for entity in participant.get("refEntities", []):
            if entity.get("schemaClass") not in {"ReferenceGeneProduct", "ReferenceIsoform"}:
                continue
            tokens = str(entity.get("displayName", "")).split()
            if len(tokens) >= 2:
                genes.add(tokens[-1])
    return genes


ANALYSIS.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)


inventory_rows = []
inventory_files = sorted(SUPPLEMENTS.glob("43856_2026_1924_MOESM*_ESM.*"))
inventory_files.append(RAW / "s43856-026-01924-9_reference.pdf")
for path in inventory_files:
    if not path.exists():
        continue
    row = {
        "file_name": path.name,
        "file_type": path.suffix.lower().lstrip("."),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "sheet_or_page_count": np.nan,
        "sheet_names": "",
        "text_extractable": np.nan,
        "role": "main_article" if path.name.endswith("reference.pdf") else "supplement",
    }
    if path.suffix.lower() == ".xlsx":
        workbook = pd.ExcelFile(path)
        row["sheet_or_page_count"] = len(workbook.sheet_names)
        row["sheet_names"] = ";".join(workbook.sheet_names)
    elif path.suffix.lower() == ".pdf":
        reader = PdfReader(path)
        row["sheet_or_page_count"] = len(reader.pages)
        row["text_extractable"] = any(bool((page.extract_text() or "").strip()) for page in reader.pages)
    inventory_rows.append(row)

inventory = pd.DataFrame(inventory_rows)
inventory.to_csv(ANALYSIS / "16_PXD083842_supplement_inventory.tsv", sep="\t", index=False)


covariate_rows = [
    {
        "domain": "demographic",
        "variable": "age_years",
        "plasma_PBC": "55 [46,61]",
        "plasma_MASLD": "49 [36,58]",
        "plasma_Healthy": "23 [21,26]",
        "plasma_group_p": "<0.001",
        "liver_PBC": "57 [51,62]",
        "liver_MASLD": "45 [34,56]",
        "liver_Healthy": "22 [21,26]",
        "liver_group_p": "<0.001",
        "patient_level_public": False,
        "TYK2_result_after_control": "not FDR-significant; no liver or plasma protein survived age adjustment",
        "risk": "critical",
        "interpretation": "Minimal age overlap makes disease and age effects non-identifiable in this small discovery cohort.",
        "source": "Main Table 1; Supplementary Table S3; Supplementary Data 3",
    },
    {
        "domain": "demographic",
        "variable": "female_sex",
        "plasma_PBC": "14/15 (93%)",
        "plasma_MASLD": "7/17 (41%)",
        "plasma_Healthy": "9/19 (47%)",
        "plasma_group_p": "0.005",
        "liver_PBC": "13/14 (93%)",
        "liver_MASLD": "6/15 (40%)",
        "liver_Healthy": "5/9 (56%)",
        "liver_group_p": "0.011",
        "patient_level_public": False,
        "TYK2_result_after_control": "not present in sex-adjusted FDR<0.05 lists",
        "risk": "high",
        "interpretation": "Sex is imbalanced and the public supplement does not expose the non-significant TYK2 coefficient.",
        "source": "Main Table 1; Supplementary Table S3; Supplementary Data 3",
    },
    {
        "domain": "metabolic",
        "variable": "BMI_kg_m2",
        "plasma_PBC": "31.7 [23.5,35.0]",
        "plasma_MASLD": "29.8 [26.6,34.7]",
        "plasma_Healthy": "22.6 [21.3,23.3]",
        "plasma_group_p": "<0.001",
        "liver_PBC": "31.8 [24.2,35.4]",
        "liver_MASLD": "29.8 [26.2,34.8]",
        "liver_Healthy": "22.6 [21.9,22.8]",
        "liver_group_p": "0.001",
        "patient_level_public": False,
        "TYK2_result_after_control": "not present in BMI-adjusted FDR<0.05 lists",
        "risk": "high",
        "interpretation": "PBC and healthy controls differ markedly in adiposity; exact adjusted TYK2 beta and P are unavailable.",
        "source": "Main Table 1; Supplementary Table S3; Supplementary Data 3",
    },
    {
        "domain": "comorbidity",
        "variable": "concurrent_MASLD_in_PBC",
        "plasma_PBC": "6/15 (40%)",
        "plasma_MASLD": "17/17",
        "plasma_Healthy": "0/19",
        "plasma_group_p": "not estimable as a covariate",
        "liver_PBC": "6/14 (43%)",
        "liver_MASLD": "15/15",
        "liver_Healthy": "0/9",
        "liver_group_p": "not estimable as a covariate",
        "patient_level_public": False,
        "TYK2_result_after_control": "not FDR-significant after excluding PBC with concurrent MASLD; not FDR-significant in PBC vs MASLD",
        "risk": "critical",
        "interpretation": "TYK2 does not currently show disease specificity relative to MASLD or after removing concurrent MASLD.",
        "source": "Main Table 1; Supplementary Table S3; Supplementary Data 2 and 4",
    },
    {
        "domain": "treatment",
        "variable": "UDCA_at_sampling",
        "plasma_PBC": "7/15 (47%)",
        "plasma_MASLD": "0/17",
        "plasma_Healthy": "0/19",
        "plasma_group_p": "<0.001",
        "liver_PBC": "6/14 (43%)",
        "liver_MASLD": "0/15",
        "liver_Healthy": "0/9",
        "liver_group_p": "0.002",
        "patient_level_public": False,
        "TYK2_result_after_control": "not evaluated for TYK2",
        "risk": "high",
        "interpretation": "Treatment is nested within PBC and cannot be separated without patient-level metadata.",
        "source": "Main Table 1; Supplementary Table S3",
    },
    {
        "domain": "disease_severity",
        "variable": "cirrhosis",
        "plasma_PBC": "1/15 (7%)",
        "plasma_MASLD": "1/17 (6%)",
        "plasma_Healthy": "0/19",
        "plasma_group_p": "0.536",
        "liver_PBC": "1/14 (7%)",
        "liver_MASLD": "1/15 (7%)",
        "liver_Healthy": "0/9",
        "liver_group_p": "0.719",
        "patient_level_public": False,
        "TYK2_result_after_control": "log2FC 0.416; FDR 0.0162 after cirrhosis adjustment",
        "risk": "low_for_TYK2",
        "interpretation": "The TYK2 association is not explained by the single cirrhotic PBC sample.",
        "source": "Main Table 1; Supplementary Table S3; Supplementary Data 3",
    },
    {
        "domain": "recruitment",
        "variable": "participant_source",
        "plasma_PBC": "FALL hospital cohort",
        "plasma_MASLD": "FALL hospital cohort",
        "plasma_Healthy": "advertisements/posters; PLS IDs",
        "plasma_group_p": "exactly confounded by public subject prefix",
        "liver_PBC": "14/14 FALL IDs",
        "liver_MASLD": "15/15 FALL IDs",
        "liver_Healthy": "9/9 PLS IDs",
        "liver_group_p": "exactly confounded by public subject prefix",
        "patient_level_public": True,
        "TYK2_result_after_control": "not controlled",
        "risk": "critical",
        "interpretation": "Disease status is inseparable from recruitment stream and subject-ID prefix in the released matrices.",
        "source": "Methods; matrix subject IDs; Reporting Summary",
    },
    {
        "domain": "biospecimen",
        "variable": "liver_biopsy_route",
        "plasma_PBC": "not applicable",
        "plasma_MASLD": "not applicable",
        "plasma_Healthy": "not applicable",
        "plasma_group_p": "not applicable",
        "liver_PBC": "percutaneous or transjugular; counts not public",
        "liver_MASLD": "percutaneous or transjugular; counts not public",
        "liver_Healthy": "9/9 liver controls underwent transjugular biopsy",
        "liver_group_p": "not calculable",
        "patient_level_public": False,
        "TYK2_result_after_control": "not evaluated",
        "risk": "high",
        "interpretation": "Biopsy route may be entangled with the healthy liver group and cannot be adjusted from public data.",
        "source": "Methods; Main Table 1 note",
    },
]

covariates = pd.DataFrame(covariate_rows)
covariates.to_csv(ANALYSIS / "17_PXD083842_covariate_metadata_audit.tsv", sep="\t", index=False)


supplement_1 = SUPPLEMENTS / "43856_2026_1924_MOESM3_ESM.xlsx"
supplement_2 = SUPPLEMENTS / "43856_2026_1924_MOESM4_ESM.xlsx"
supplement_3 = SUPPLEMENTS / "43856_2026_1924_MOESM5_ESM.xlsx"
supplement_4 = SUPPLEMENTS / "43856_2026_1924_MOESM6_ESM.xlsx"

analysis_definitions = [
    {
        "analysis": "primary_unadjusted",
        "comparison": "PBC_vs_Healthy",
        "adjustment_or_exclusion": "none",
        "path": supplement_1,
        "up": "upreg_liver_pbc_vs_healthy",
        "down": "downreg_liver_pbc_vs_healthy",
    },
    {
        "analysis": "exclude_concurrent_MASLD",
        "comparison": "PBC_without_MASLD_vs_Healthy",
        "adjustment_or_exclusion": "exclude PBC with concurrent MASLD",
        "path": supplement_2,
        "up": "upreg_liver_pbc_vs_healthy",
        "down": "downreg_liver_pbc_vs_healthy",
    },
    {
        "analysis": "age_adjusted",
        "comparison": "PBC_vs_Healthy",
        "adjustment_or_exclusion": "age",
        "path": supplement_3,
        "up": "liver_age_up",
        "down": "liver_age_down",
    },
    {
        "analysis": "sex_adjusted",
        "comparison": "PBC_vs_Healthy",
        "adjustment_or_exclusion": "sex",
        "path": supplement_3,
        "up": "liver_sex_up",
        "down": "liver_sex_down",
    },
    {
        "analysis": "BMI_adjusted",
        "comparison": "PBC_vs_Healthy",
        "adjustment_or_exclusion": "BMI",
        "path": supplement_3,
        "up": "liver_bmi_up",
        "down": "liver_bmi_down",
    },
    {
        "analysis": "cirrhosis_adjusted",
        "comparison": "PBC_vs_Healthy",
        "adjustment_or_exclusion": "cirrhosis",
        "path": supplement_3,
        "up": "liver_cirrhosis_up",
        "down": "liver_cirrhosis_down",
    },
    {
        "analysis": "age_plus_sex_adjusted",
        "comparison": "PBC_vs_Healthy",
        "adjustment_or_exclusion": "age + sex",
        "path": supplement_3,
        "up": "liver_age_sex_up",
        "down": "liver_age_sex_down",
    },
    {
        "analysis": "disease_specific_comparator",
        "comparison": "PBC_vs_MASLD",
        "adjustment_or_exclusion": "direct disease comparator",
        "path": supplement_4,
        "up": "upreg_liver_pbc_vs_masld",
        "down": "downreg_liver_pbc_vs_masld",
    },
]

tyk2_rows = []
ifng_rows = []
liver, liver_groups, liver_subjects = read_matrix(RAW / "liver_protein_data_minprob_processed.tsv.txt")
ifng_genes = reactome_genes(RAW / "reactome_participants_R-HSA-877300.json")
ifng_measured = sorted(gene for gene in ifng_genes if gene in liver.index and not gene.startswith("HLA-") and gene != "B2M")

for definition in analysis_definitions:
    up_genes, down_genes = sheet_gene_sets(definition["path"], definition["up"], definition["down"])
    up_row_count = len(significant_sheet(definition["path"], definition["up"]))
    down_row_count = len(significant_sheet(definition["path"], definition["down"]))
    tyk2_row, tyk2_sheet = first_gene_row(definition["path"], [definition["up"], definition["down"]], "TYK2")
    exact = tyk2_row is not None
    tyk2_rows.append(
        {
            "dataset": "PXD083842",
            "tissue": "liver",
            "analysis": definition["analysis"],
            "comparison": definition["comparison"],
            "adjustment_or_exclusion": definition["adjustment_or_exclusion"],
            "significant_protein_rows_up_n": up_row_count,
            "significant_protein_rows_down_n": down_row_count,
            "TYK2_in_FDR_lt_0_05_list": exact,
            "direction": "up" if exact and tyk2_sheet == definition["up"] else ("down" if exact else "not_reported"),
            "log2FC": float(tyk2_row["log2FC"]) if exact else np.nan,
            "adjusted_p": float(tyk2_row["adj.P.Val"]) if exact else np.nan,
            "exact_non_significant_beta_p_public": False if not exact else True,
            "interpretation": (
                "FDR-significant in the authors' released table"
                if exact
                else "Absent from the authors' FDR<0.05-only table; exact beta and P are not public"
            ),
        }
    )
    significant_ifng_up = sorted(set(ifng_measured) & up_genes)
    significant_ifng_down = sorted(set(ifng_measured) & down_genes)
    ifng_rows.append(
        {
            "dataset": "PXD083842",
            "tissue": "liver",
            "analysis": definition["analysis"],
            "comparison": definition["comparison"],
            "adjustment_or_exclusion": definition["adjustment_or_exclusion"],
            "prespecified_IFNG_non_HLA_B2M_measured_n": len(ifng_measured),
            "significant_up_n": len(significant_ifng_up),
            "significant_down_n": len(significant_ifng_down),
            "significant_up_genes": ";".join(significant_ifng_up),
            "significant_down_genes": ";".join(significant_ifng_down),
            "can_recompute_adjusted_pathway_score": False,
            "limitation": "Only FDR<0.05 protein lists are public; full adjusted coefficients are unavailable.",
        }
    )

tyk2_results = pd.DataFrame(tyk2_rows)
tyk2_results.to_csv(ANALYSIS / "18_PXD083842_TYK2_adjusted_results.tsv", sep="\t", index=False)
pd.DataFrame(ifng_rows).to_csv(ANALYSIS / "18b_PXD083842_IFNG_adjusted_significant_members.tsv", sep="\t", index=False)


plasma_groups, plasma_subjects = read_matrix_metadata(RAW / "plasma_protein_data_minprob_processed.tsv.txt")
all_subjects = sorted(set(liver_subjects) | set(plasma_subjects))
paired_subjects = sorted(set(liver_subjects) & set(plasma_subjects))
plasma_only_subjects = sorted(set(plasma_subjects) - set(liver_subjects))

sample_number = pd.Series(
    [int(re.search(r"\d+", sample).group()) for sample in liver_groups.index],
    index=liver_groups.index,
    name="sample_number_proxy",
)
tyk2 = liver.loc["TYK2"].astype(float)
sample_number_test = stats.kruskal(*[sample_number[liver_groups == group] for group in ["PBC", "NAFLD", "Healthy"]])
sample_number_tyk2 = stats.spearmanr(sample_number, tyk2)

batch_rows = [
    {
        "audit_item": "sample_tissue_counts",
        "value": "liver 38 (14 PBC, 15 MASLD/NAFLD, 9 healthy); plasma 51 (15 PBC, 17 MASLD/NAFLD, 19 healthy)",
        "status": "verified_from_matrices",
        "risk": "low",
        "interpretation": "Counts agree with the published article after accounting for one PBC, two MASLD and ten healthy plasma-only participants.",
    },
    {
        "audit_item": "pairing",
        "value": f"{len(paired_subjects)}/38 liver participants have matched plasma; {len(plasma_only_subjects)}/51 plasma participants lack liver",
        "status": "verified_from_subject_IDs",
        "risk": "moderate_for_cross_tissue_claims",
        "interpretation": "All released liver samples are paired to plasma, but the full plasma cohort is not fully paired.",
    },
    {
        "audit_item": "participant_source_prefix",
        "value": "all PBC/MASLD are FALL IDs; all healthy controls are PLS IDs in both tissues",
        "status": "exact_confounding",
        "risk": "critical",
        "interpretation": "Recruitment stream/source cannot be separated from disease group in the released data.",
    },
    {
        "audit_item": "analytical_batch_statement",
        "value": "no batch correction; all samples from each tissue were analysed together in a single analytical run",
        "status": "reported_by_authors",
        "risk": "moderate",
        "interpretation": "This reduces discrete batch concerns but does not exclude within-run order drift or preparation-day effects.",
    },
    {
        "audit_item": "run_randomization_and_QC",
        "value": "injection order, sample randomization, pooled-QC placement and drift model not reported",
        "status": "not_public",
        "risk": "high",
        "interpretation": "The reporting summary states participants were not randomized; it does not document analytical run-order randomization.",
    },
    {
        "audit_item": "blinding",
        "value": "not performed",
        "status": "reported_by_authors",
        "risk": "moderate",
        "interpretation": "Blinding is less central for automated DIA acquisition but remains relevant for sample handling and exclusions.",
    },
    {
        "audit_item": "raw_file_completeness",
        "value": "1 raw DIA archive for 89 sample-tissue measurements (1.12% if one raw file per measurement)",
        "status": "PRIDE_PARTIAL_submission",
        "risk": "critical_for_batch_reanalysis",
        "interpretation": "Full run-order, QC, missingness and batch reprocessing cannot be reconstructed from PRIDE.",
    },
    {
        "audit_item": "sample_number_proxy_by_group",
        "value": f"liver Kruskal P={sample_number_test.pvalue:.6g}; healthy S54-S62, PBC S01-S50, MASLD S04-S47",
        "status": "strong_group_separation_but_not_confirmed_run_order",
        "risk": "high_if_sample_ID_reflects_acquisition_order",
        "interpretation": "Treat as a warning only: sample numbers are not explicitly defined as injection order.",
    },
    {
        "audit_item": "TYK2_vs_sample_number_proxy",
        "value": f"Spearman rho={sample_number_tyk2.statistic:.3f}, P={sample_number_tyk2.pvalue:.4g}",
        "status": "association_driven_by_between_group_structure",
        "risk": "warning_not_proof_of_drift",
        "interpretation": "Within-group correlations are non-significant; the overall correlation cannot distinguish biology from group/order confounding.",
    },
    {
        "audit_item": "processed_matrix",
        "value": "global median normalization; log2; >=70% valid in at least one group; MinProb width 0.3/downshift 1.8",
        "status": "reported_by_authors",
        "risk": "high_for_missingness_sensitivity",
        "interpretation": "Unimputed intensities and full raw files are absent, so the TYK2 missingness mechanism cannot be audited.",
    },
]

batch_pairing = pd.DataFrame(batch_rows)
batch_pairing.to_csv(ANALYSIS / "19_PXD083842_batch_pairing_audit.tsv", sep="\t", index=False)


primary = tyk2_results.set_index("analysis").loc["primary_unadjusted"]
cirrhosis = tyk2_results.set_index("analysis").loc["cirrhosis_adjusted"]
report = f"""# PXD083842 样本元数据、批次与混杂因素审计

## 审计结论

`PXD083842` 已有正式同行评议论文：Kellemann 等，*Communications Medicine*，2026-10-01，DOI `10.1038/s43856-026-01924-9`。它不再属于“无关联论文”的未发表资源，但仍是 `PARTIAL` PRIDE submission。

肝组织 TYK2 的未校正差异确实存在，作者原始 limma 结果为 log2FC={primary['log2FC']:.3f}、FDR={primary['adjusted_p']:.4g}；肝硬化校正后仍为 log2FC={cirrhosis['log2FC']:.3f}、FDR={cirrhosis['adjusted_p']:.4g}。这与此前的逐样本删除、置换和稳健回归结果一致，说明结果不是由单一离群样本或唯一肝硬化病例制造。

但是，TYK2 **未进入**年龄校正、年龄+性别校正、性别校正、BMI校正、排除合并MASLD或PBC-vs-MASLD的FDR<0.05列表。尤其是年龄校正后，作者报告肝和血浆均无任何蛋白通过FDR。由于健康对照极年轻，年龄在肝组织中为 PBC 57 [51,62] 岁、healthy 22 [21,26] 岁，组间重叠极小，当前公开数据不能把年龄效应与PBC效应可靠拆开。补充表只发布显著蛋白，未发布TYK2未显著模型的beta、SE和P，因此不能把“未进入列表”解释为方向反转。

预设的非HLA/B2M IFN-γ程序比单个TYK2更广泛：48个实测成员中，未校正时有19个上调蛋白通过FDR，排除合并MASLD后仍有8个，性别校正12个、BMI校正10个、肝硬化校正18个，PBC-vs-MASLD仍有10个；年龄或年龄+性别校正后为0。由于作者只公开FDR显著列表，不能复算完整校正后的通路分数。最稳妥的解释是“免疫程序在多种非年龄敏感性分析中仍可见，但年龄与疾病仍难分离”，而不是把通路显著当作TYK2本身的校正复制。

**证据等级调整：** PXD083842 对 TYK2 应从“独立强蛋白复制”降为“方向一致、离群值稳健，但临床混杂明显的支持性肝蛋白证据”。它不能单独证明PBC特异性，也不能作为TYK2药物再利用的决定性证据。

## 主要混杂

- 年龄、性别和BMI均显著失衡；年龄是最严重且近乎不可辨识的混杂。
- PBC中40%（全血浆队列）/43%（肝子集）合并MASLD；排除这些患者后TYK2不再进入FDR显著列表。
- PBC中47%（全队列）/43%（肝子集）采样时已使用UDCA，公开患者级数据不足以校正。
- 所有PBC/MASLD受试者均为 `FALL_*`，所有健康对照均为 `PLS_*`；招募来源与疾病分组完全重合。
- 9名具有肝组织的健康对照均经颈静脉肝活检；患者的活检路径分布未公开，存在组织获取方式混杂。

## 批次与原始数据

作者说明每种组织的样本在一个连续分析批次中完成，因此未做批次校正。公开方法未报告注射顺序随机化、 pooled QC位置、漂移曲线或样本制备批次。肝样本编号在组间明显分离（编号代理的Kruskal P={sample_number_test.pvalue:.3g}），且TYK2与编号代理相关（rho={sample_number_tyk2.statistic:.3f}，P={sample_number_tyk2.pvalue:.4g}），但编号并未被作者定义为注射顺序，所以这只是风险警报，不能当作已证实批次效应。

PRIDE只公开1个原始DIA压缩包，而两个矩阵合计有89个样本-组织测量；因此无法从公开资源重建完整运行顺序、QC、缺失值分布或进行无填补复算。MinProb填补后的矩阵仍可做外部方向核验，但不适合宣称完成了完整技术批次校正。

## 配对结构

- 肝组织38人：14 PBC、15 MASLD、9 healthy。
- 血浆51人：15 PBC、17 MASLD、19 healthy。
- 全部38名肝组织受试者都有对应血浆；另有13名仅血浆受试者（1 PBC、2 MASLD、10 healthy）。
- 正确表述应为：“肝组织子集均有配对血浆，但完整血浆队列并非人人有肝组织”。血浆中未定量TYK2，因此不能声称TYK2跨组织复制。

## 对主线的影响

1. 保留PXD083842，但把它放在遗传方向、独立GWAS和外部转录组证据之后。
2. 主文措辞限定为“PBC相对年轻健康肝组织中TYK2蛋白升高；该关联对离群值稳健，但在年龄/代谢敏感性分析中减弱”。
3. 不写“PBC特异蛋白标志物”，不写“经完整多变量校正仍显著”，不把肝TYK2写成血浆复制。
4. 若要真正解决混杂，需要作者提供患者级年龄、性别、BMI、MASLD、UDCA、活检路径、注射顺序、QC和未填补强度矩阵；现有公开数据不能通过放宽阈值补救这一识别问题。

## 输出

- `16_PXD083842_supplement_inventory.tsv`
- `17_PXD083842_covariate_metadata_audit.tsv`
- `18_PXD083842_TYK2_adjusted_results.tsv`
- `18b_PXD083842_IFNG_adjusted_significant_members.tsv`
- `19_PXD083842_batch_pairing_audit.tsv`
"""

(REPORTS / "03_PXD083842_confounder_audit_report.md").write_text(report, encoding="utf-8")

print(f"Inventory files: {len(inventory)}")
print(f"Covariates audited: {len(covariates)}")
print(f"TYK2 analyses audited: {len(tyk2_results)}")
print(f"Paired liver-plasma subjects: {len(paired_subjects)}")
print(f"Plasma-only subjects: {len(plasma_only_subjects)}")
