from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Polygon, Patch
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "figures" / "main"
ANALYSIS = REPO / "data" / "derived" / "analysis" / "single_cell"
CARDON = ANALYSIS / "cardon_antigen_context"
EXTERNAL = ANALYSIS / "external"
AUDIT = ANALYSIS / "resource_audit" / "independent_pbc_dataset_audit.tsv"

ANTIGEN = CARDON / "02_pure_PBC_patient_condition_pseudobulk.tsv"
ANTIGEN_SUMMARY = CARDON / "04_paired_inference_summary.tsv"
COHORT = EXTERNAL / "04_cohort_specificity_tests.tsv"
META = EXTERNAL / "05_cross_cohort_specificity_meta.tsv"
DECISION = EXTERNAL / "07_external_scRNA_replication_decision.tsv"

INK = "#263746"
MUTED = "#718096"
GRID = "#DCE4EA"
BLUE = "#2B7BBB"
SKY = "#70B7E6"
GREEN = "#2A9D8F"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
GRAY = "#AAB5BC"
PALE = "#EEF3F6"


def read(path):
    return pd.read_csv(path, sep="\t")


def style_axis(ax, grid_axis=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.7, alpha=0.9)
    ax.set_axisbelow(True)


def panel_label(ax, label, x=-0.09, y=1.08):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=14, fontweight="bold", color=INK, va="top")


def audit_matrix():
    audit = read(AUDIT).set_index("accession")
    rows = [
        "HRA008003",
        "HRA002347",
        "CELLxGENE_J_Hepatol_2024",
        "GSE243981",
        "GSE270739",
        "GSE281263",
        "JHEP_2025_CCR5_Trm",
        "Gut_2025_Kupffer_autophagy",
    ]
    labels = [
        "HRA008003",
        "HRA002347",
        "J Hepatol atlas",
        "GSE243981",
        "GSE270739",
        "GSE281263",
        "CCR5 Trm study",
        "Kupffer study",
    ]
    # 2=eligible/yes, 1=limited/context/unclear, 0=no.
    values = np.array([
        [2, 2, 2, 2, 2, 2],  # HRA008003, primary cohort
        [2, 2, 2, 2, 2, 1],  # paired platforms; external liver context
        [2, 2, 2, 2, 2, 2],  # donor-level public atlas
        [2, 2, 2, 0, 2, 0],  # overlaps the atlas framework
        [0, 2, 2, 2, 2, 0],  # valid paired antigen context, not case-control
        [0, 2, 1, 2, 1, 0],  # technical/condition units, inadequate PBC donor replication
        [2, 1, 0, 1, 1, 0],  # relevant tissue, no reusable public matrix
        [2, 1, 0, 1, 1, 0],  # relevant tissue, no reusable accession
    ], dtype=int)
    columns = ["PBC liver", "Donor IDs", "Reusable\nmatrix", "Independent\ndonors", "IL12RB2\nmeasurable", "Valid PBC–control\ncontrast"]
    decisions = [audit.loc[r, "decision"] for r in rows]
    return labels, columns, values, decisions


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.4,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.0,
        "legend.fontsize": 7.2,
        "legend.frameon": False,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })

    antigen = read(ANTIGEN)
    antigen_summary = read(ANTIGEN_SUMMARY)
    cohort = read(COHORT)
    meta = read(META).iloc[0]

    fig = plt.figure(figsize=(14.7, 8.5), facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.03, 1.18], height_ratios=[1.0, 1.02], wspace=0.32, hspace=0.44)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    # A: paired, patient-level antigen context. This is deliberately not a UMAP.
    p = antigen.pivot(index="patient_core", columns="condition", values="IL12RB2_log2CPM").sort_index()
    x0, x1 = 0, 1
    endpoint_rows = []
    for patient, row in p.iterrows():
        y0 = row["foreign antigen"]
        y1 = row["PDC-E2"]
        color = GREEN if y1 > y0 else GRAY
        ax_a.plot([x0, x1], [y0, y1], color=color, alpha=0.72, linewidth=1.6, zorder=1)
        ax_a.scatter([x0, x1], [y0, y1], color=[BLUE, VERMILLION], s=38, edgecolor="white", linewidth=0.6, zorder=2)
        endpoint_rows.append((patient, y1))
    # Deterministic label repulsion for the three close high-expression endpoints.
    endpoint_rows.sort(key=lambda x: x[1], reverse=True)
    label_positions = []
    for _, observed_y in endpoint_rows:
        label_y = observed_y if not label_positions else min(observed_y, label_positions[-1] - 0.37)
        label_positions.append(label_y)
    for (patient, observed_y), label_y in zip(endpoint_rows, label_positions):
        ax_a.annotate(patient, xy=(1.0, observed_y), xytext=(1.08, label_y),
                      va="center", fontsize=6.5, color=MUTED,
                      arrowprops=dict(arrowstyle="-", color="#BAC3C9", lw=0.65))
    summary = antigen_summary[(antigen_summary.metric == "IL12RB2_log2CPM") & (antigen_summary.analysis == "pure_PBC")].iloc[0]
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["Foreign antigen", "PDC-E2"])
    ax_a.set_xlim(-0.22, 1.43)
    ax_a.set_ylabel("IL12RB2 patient pseudobulk\nlog2(CPM + 0.5)")
    ax_a.set_title("Within-PBC antigen-context test", loc="left", color=INK)
    ax_a.text(0.02, 0.97,
              f"5 paired PBC patients\nMean difference={summary.mean_paired_difference:.3f}\n95% CI {summary.ci95_low:.3f} to {summary.ci95_high:.3f}\nExact sign-flip P={summary.exact_signflip_p_two_sided:.3f}",
              transform=ax_a.transAxes, va="top", fontsize=7.2, color=INK,
              bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=GRID, alpha=0.94))
    ax_a.text(0.99, 0.02, "Antigen context only; not PBC vs healthy", transform=ax_a.transAxes, ha="right", fontsize=6.8, color=VERMILLION)
    style_axis(ax_a, "y")
    panel_label(ax_a, "A")

    # B: cohort-level lrNK-cNK estimates. Three independent cohorts only.
    b = cohort[(cohort.outcome == "logCPM") & (cohort.independent_primary == True)].copy()
    order = ["CELLxGENE_scRNA", "HRA002347_5prime", "HRA008003"]
    b["sort"] = b.cohort.map({v: i for i, v in enumerate(order)})
    b = b.sort_values("sort")
    names = {
        "CELLxGENE_scRNA": "J Hepatol atlas",
        "HRA002347_5prime": "HRA002347 5′",
        "HRA008003": "HRA008003",
    }
    colors = {"CELLxGENE_scRNA": BLUE, "HRA002347_5prime": ORANGE, "HRA008003": GREEN}
    ypos = np.arange(len(b))[::-1]
    for y, (_, row) in zip(ypos, b.iterrows()):
        color = colors[row.cohort]
        ax_b.errorbar(row.mean_contrast, y,
                      xerr=[[row.mean_contrast - row.ci_low], [row.ci_high - row.mean_contrast]],
                      fmt="o", color=color, markersize=6.4, capsize=3.0, linewidth=1.5)
        ax_b.text(3.52, y, f"n={int(row.donors)}; sign-flip P={row.exact_sign_flip_p_one_sided:.4f}", va="center", fontsize=6.8, color=INK)
    ax_b.axvline(0, color="#66717A", linewidth=0.9)
    ax_b.set_yticks(ypos)
    ax_b.set_yticklabels([names[x] for x in b.cohort])
    ax_b.set_xlim(-0.65, 6.0)
    ax_b.set_xlabel("Within-donor IL12RB2 logCPM: lrNK − cNK")
    ax_b.set_title("Non-overlapping liver cellular-context cohorts", loc="left", color=INK)
    ax_b.text(0.02, 0.04, "Non-overlapping primary cohorts only", transform=ax_b.transAxes,
              ha="left", fontsize=6.8, color=MUTED,
              bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.88))
    style_axis(ax_b, "x")
    panel_label(ax_b, "B")

    # C: the random-effects conclusion as a separate evidence panel.
    y = 0.64
    ax_c.errorbar(meta.fixed_effect, y + 0.17,
                  xerr=[[meta.fixed_effect - meta.fixed_effect_ci_low], [meta.fixed_effect_ci_high - meta.fixed_effect]],
                  fmt="o", color=GRAY, markersize=5.5, capsize=2.5, label="Fixed effect")
    diamond = Polygon([
        [meta.random_effect_HK_ci_low, y - 0.17],
        [meta.random_effect_DL, y - 0.07],
        [meta.random_effect_HK_ci_high, y - 0.17],
        [meta.random_effect_DL, y - 0.27],
    ], closed=True, facecolor=ORANGE, edgecolor="#A76F00", linewidth=0.9)
    ax_c.add_patch(diamond)
    ax_c.axvline(0, color="#66717A", linewidth=0.9)
    ax_c.set_yticks([y + 0.17, y - 0.17])
    ax_c.set_yticklabels(["Fixed effect", "HK random effects"])
    ax_c.set_ylim(0.08, 1.02)
    ax_c.set_xlim(-2.2, 5.7)
    ax_c.set_xlabel("Pooled lrNK − cNK IL12RB2 logCPM")
    ax_c.set_title("Hartung–Knapp cross-cohort synthesis", loc="left", color=INK)
    ax_c.text(0.03, 0.08,
              f"3 cohorts · {int(meta.donors)} donors · 3/3 positive\nHK estimate={meta.random_effect_DL:.3f}\n95% CI {meta.random_effect_HK_ci_low:.3f} to {meta.random_effect_HK_ci_high:.3f}\nOne-sided P={meta.random_effect_HK_p_one_sided:.3f}\nI²={meta.I2_percent:.1f}% · τ²={meta.tau2_DL:.3f}",
              transform=ax_c.transAxes, va="bottom", fontsize=7.2, color=INK,
              bbox=dict(boxstyle="round,pad=0.32", facecolor="white", edgecolor=GRID, alpha=0.95))
    ax_c.text(0.98, 0.08, "Directionally consistent\nbut highly heterogeneous", transform=ax_c.transAxes,
              ha="right", va="bottom", fontsize=8.0, color=VERMILLION, fontweight="bold")
    style_axis(ax_c, "x")
    panel_label(ax_c, "C")

    # D: systematic resource audit; no fourth eligible donor-level case-control cohort.
    labels, columns, values, decisions = audit_matrix()
    cmap = ListedColormap(["#D6DDE2", "#F1C76A", "#3BAF8B"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    ax_d.imshow(values, cmap=cmap, norm=norm, aspect="auto")
    ax_d.set_xticks(np.arange(len(columns)))
    ax_d.set_xticklabels(columns, rotation=0, ha="center")
    ax_d.set_yticks(np.arange(len(labels)))
    ax_d.set_yticklabels(labels)
    ax_d.tick_params(length=0)
    for i in range(values.shape[0] + 1):
        ax_d.axhline(i - 0.5, color="white", linewidth=1.4)
    for j in range(values.shape[1] + 1):
        ax_d.axvline(j - 0.5, color="white", linewidth=1.4)
    symbols = {0: "×", 1: "△", 2: "✓"}
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax_d.text(j, i, symbols[values[i, j]], ha="center", va="center",
                      fontsize=9.5, fontweight="bold", color="white" if values[i, j] == 2 else INK)
    ax_d.set_title("Systematic audit of public PBC single-cell resources", loc="left", color=INK)
    legend = [Patch(facecolor="#3BAF8B", label="Yes/eligible"), Patch(facecolor="#F1C76A", label="Limited/context only"), Patch(facecolor="#D6DDE2", label="No/unavailable")]
    ax_d.legend(handles=legend, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.14), columnspacing=1.1, handlelength=1.0)
    ax_d.text(0.99, -0.29, "Audit updated: Oct 20, 2026 · No additional adequately replicated donor-level PBC–control cohort identified",
              transform=ax_d.transAxes, ha="right", fontsize=6.8, color=VERMILLION, fontweight="bold")
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    panel_label(ax_d, "D")

    fig.suptitle("Figure 6 | External single-cell resource audit and cellular-context synthesis",
                 x=0.045, y=0.992, ha="left", fontsize=15.6, fontweight="bold", color=INK)
    fig.text(0.995, 0.005,
             "Donor-level evidence supports recurrent lrNK enrichment of IL12RB2; small cohort numbers and high heterogeneity preclude definitive replication.",
             ha="right", fontsize=7.0, color=MUTED)
    fig.savefig(OUT / "Fig06_external_single_cell_context_locked.pdf")
    fig.savefig(OUT / "Fig06_external_single_cell_context_locked.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
