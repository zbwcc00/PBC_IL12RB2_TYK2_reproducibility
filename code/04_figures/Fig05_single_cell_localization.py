from pathlib import Path
import gzip
import math

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import norm


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "figures" / "main"
INPUT = REPO / "data" / "derived" / "figures" / "main" / "Fig05"
HRA = REPO / "data" / "derived" / "analysis" / "single_cell" / "hra008003" / "pseudobulk"

UMAP = INPUT / "01_HRA008003_UMAP_and_targets_deidentified.tsv.gz"
MARKERS = INPUT / "02_marker_dotplot_summary.tsv"
TARGET_DOT = INPUT / "03_target_axis_dotplot_summary.tsv"
COMPOSITION = INPUT / "04_donor_NK_composition.tsv"
COMPOSITION_TEST = INPUT / "05_NK_composition_exact_permutation.tsv"
EDGE = HRA / "03_prespecified_target_edgeR_results.tsv"
PERM = HRA / "05_exact_donor_permutation_results.tsv"
LODO = HRA / "07_target_leave_one_donor_out.tsv"

INK = "#263746"
MUTED = "#718096"
GRID = "#DCE4EA"
BLUE = "#2B7BBB"
SKY = "#70B7E6"
TEAL = "#16A085"
GREEN = "#2A9D8F"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
PURPLE = "#8C6BB1"
PINK = "#CC79A7"
LIGHT = "#ECF2F6"

CELL_COLORS = {
    "Hepatocyte": "#D8A43A",
    "Cholangiocyte": "#E76F51",
    "Endothelial": "#5AAE9C",
    "Fibroblast_Stellate": "#B58AC4",
    "Myeloid": "#7A8D9B",
    "Neutrophil": "#A8B3BA",
    "B_cell": "#4C78A8",
    "Plasma_cell": "#8CB7D9",
    "CD4_T": "#009E73",
    "CD8_T": "#0072B2",
    "lrNK": "#D55E00",
    "cNK": "#E69F00",
    "T_other": "#6C8E67",
    "Proliferating": "#CC79A7",
    "RBC": "#B23A48",
    "Mast_cell": "#9C755F",
    "uncertain": "#C6CDD2",
    "TNK_mixed": "#536D7A",
}


def style_axis(ax, grid_axis=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.65, alpha=0.9)
    ax.set_axisbelow(True)


def panel_label(ax, label, x=-0.10, y=1.08):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=14, fontweight="bold", color=INK, va="top")


def read(path):
    return pd.read_csv(path, sep="\t")


def ci_from_f(beta, f_value):
    se = abs(beta) / math.sqrt(f_value) if f_value > 0 else np.nan
    return beta - 1.96 * se, beta + 1.96 * se


def scale_within_gene(df, column):
    result = df.copy()
    result["scaled"] = result.groupby("gene")[column].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0
    )
    return result


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.2,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.7,
        "legend.fontsize": 7.0,
        "legend.frameon": False,
        "xtick.labelsize": 7.2,
        "ytick.labelsize": 7.2,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })

    umap = read(UMAP)
    marker = read(MARKERS)
    target_dot = read(TARGET_DOT)
    composition = read(COMPOSITION)
    composition_test = read(COMPOSITION_TEST)
    edge = read(EDGE)
    permutation = read(PERM)
    lodo = read(LODO)

    fig = plt.figure(figsize=(15.2, 10.0), facecolor="white")
    outer = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.05, 1.12], height_ratios=[1.02, 1.0], wspace=0.43, hspace=0.42)
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    c_grid = outer[0, 2].subgridspec(2, 1, height_ratios=[1.20, 0.82], hspace=0.34)
    ax_c1 = fig.add_subplot(c_grid[0, 0])
    ax_c2 = fig.add_subplot(c_grid[1, 0])
    ax_d = fig.add_subplot(outer[1, 0])
    ax_e = fig.add_subplot(outer[1, 1])
    ax_f = fig.add_subplot(outer[1, 2])

    # A: one donor-aware UMAP, with only large groups directly labeled.
    rng = np.random.default_rng(51005)
    order = rng.permutation(len(umap))
    for label, group in umap.iloc[order].groupby("display_label", sort=False):
        ax_a.scatter(group.UMAP1, group.UMAP2, s=1.05, alpha=0.55, linewidths=0, rasterized=True,
                     c=CELL_COLORS.get(label, "#AAB4BB"), label=label)
    label_groups = ["Hepatocyte", "Cholangiocyte", "Endothelial", "Fibroblast_Stellate", "Myeloid", "B_cell", "CD4_T", "CD8_T", "lrNK", "cNK"]
    label_text = {"Fibroblast_Stellate": "Stellate", "Cholangiocyte": "Cholangio.", "Endothelial": "Endothelial", "Hepatocyte": "Hepatocyte", "Myeloid": "Myeloid", "B_cell": "B", "CD4_T": "CD4 T", "CD8_T": "CD8 T", "lrNK": "lrNK", "cNK": "cNK"}
    label_offset = {"Hepatocyte": (-0.45, 0.65), "Cholangiocyte": (0.55, -0.55), "Fibroblast_Stellate": (0.55, 0.15), "cNK": (0.0, 0.65)}
    for label in label_groups:
        group = umap[umap.display_label == label]
        if len(group) == 0:
            continue
        x, y = group[["UMAP1", "UMAP2"]].median()
        dx, dy = label_offset.get(label, (0, 0))
        x, y = x + dx, y + dy
        ax_a.text(x, y, label_text[label], ha="center", va="center", fontsize=6.2, fontweight="bold", color=INK,
                  bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.72))
    ax_a.set_xlabel("UMAP 1")
    ax_a.set_ylabel("UMAP 2")
    ax_a.set_title("Donor-aware liver atlas", loc="left", color=INK)
    ax_a.text(0.01, 0.01, "69,436 cells · 5 PBC + 5 control donors", transform=ax_a.transAxes, fontsize=7.0, color=MUTED)
    style_axis(ax_a)
    panel_label(ax_a, "A")

    # B: compact marker-supported annotation dot plot.
    selected_markers = ["ALB", "KRT19", "PECAM1", "COL1A1", "LST1", "MS4A1", "MZB1", "CD3D", "CD8A", "CXCR6", "FCGR3A"]
    cell_order = ["Hepatocyte", "Cholangiocyte", "Endothelial", "Fibroblast_Stellate", "Myeloid", "B_cell", "Plasma_cell", "CD4_T", "CD8_T", "lrNK", "cNK"]
    b = marker[marker.gene.isin(selected_markers) & marker.cell_type.isin(cell_order)].copy()
    b = scale_within_gene(b, "mean_logNorm")
    for _, row in b.iterrows():
        x = selected_markers.index(row.gene)
        y = cell_order.index(row.cell_type)
        ax_b.scatter(x, y, s=8 + row.pct_positive * 4.0, c=row.scaled, cmap="RdYlBu_r", vmin=-1.8, vmax=1.8,
                     edgecolor=INK if row.expected_marker else "white", linewidth=0.35 if row.expected_marker else 0.15)
    ax_b.set_xticks(range(len(selected_markers)))
    ax_b.set_xticklabels(selected_markers, rotation=55, ha="right", fontstyle="italic")
    ax_b.set_yticks(range(len(cell_order)))
    ax_b.set_yticklabels(["Hepatocyte", "Cholangiocyte", "Endothelial", "Stellate", "Myeloid", "B cell", "Plasma", "CD4 T", "CD8 T", "lrNK", "cNK"])
    ax_b.invert_yaxis()
    ax_b.set_xlim(-0.7, len(selected_markers) - 0.3)
    ax_b.set_title("SingleR + canonical-marker annotation", loc="left", color=INK)
    ax_b.text(0.01, -0.24, "Size: expressing cells    color: within-gene mean expression", transform=ax_b.transAxes, fontsize=6.7, color=MUTED)
    ax_b.tick_params(length=0)
    for x in range(len(selected_markers) + 1):
        ax_b.axvline(x - 0.5, color="#EEF2F5", lw=0.45, zorder=0)
    for y in range(len(cell_order) + 1):
        ax_b.axhline(y - 0.5, color="#EEF2F5", lw=0.45, zorder=0)
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    panel_label(ax_b, "B")

    # C: one feature UMAP + four-gene T/NK dot plot.
    cmap = LinearSegmentedColormap.from_list("il12", ["#E6EBEF", "#7EC7C4", "#E76F51"])
    nonzero = umap[umap.IL12RB2_logNorm > 0].copy()
    background = umap[umap.IL12RB2_logNorm == 0]
    ax_c1.scatter(background.UMAP1, background.UMAP2, s=0.65, c="#DDE3E7", alpha=0.34, linewidths=0, rasterized=True)
    sc = ax_c1.scatter(nonzero.UMAP1, nonzero.UMAP2, s=3.4, c=nonzero.IL12RB2_logNorm, cmap=cmap,
                      vmin=0, vmax=max(1.2, nonzero.IL12RB2_logNorm.quantile(0.98)), linewidths=0, alpha=0.9, rasterized=True)
    ax_c1.set_title("IL12RB2 and pathway localization", loc="left", color=INK)
    ax_c1.set_xlabel("UMAP 1")
    ax_c1.set_ylabel("UMAP 2")
    ax_c1.text(0.99, 0.02, f"{100*(umap.IL12RB2_UMI>0).mean():.1f}% positive", transform=ax_c1.transAxes, ha="right", fontsize=6.8, color=MUTED)
    style_axis(ax_c1)
    panel_label(ax_c1, "C")

    c2 = scale_within_gene(target_dot, "mean_logNorm")
    genes = ["IL12RB2", "TYK2", "STAT4", "IFNG"]
    celltypes = ["lrNK", "cNK", "CD4_T", "CD8_T"]
    for _, row in c2.iterrows():
        x = genes.index(row.gene)
        y = celltypes.index(row.cell_type)
        ax_c2.scatter(x, y, s=18 + row.pct_positive * 4.0, c=row.scaled, cmap="PuOr_r", vmin=-1.7, vmax=1.7,
                      edgecolor="white", linewidth=0.4)
    ax_c2.set_xticks(range(4))
    ax_c2.set_xticklabels(genes, rotation=25, ha="right", fontstyle="italic")
    ax_c2.set_yticks(range(4))
    ax_c2.set_yticklabels(["lrNK", "cNK", "CD4 T", "CD8 T"])
    ax_c2.invert_yaxis()
    ax_c2.tick_params(length=0)
    ax_c2.text(1.0, -0.31, "Size: % positive · color: within-gene mean", transform=ax_c2.transAxes, ha="right", fontsize=6.3, color=MUTED)
    for spine in ax_c2.spines.values():
        spine.set_visible(False)

    # D: donor-level composition (no cell-level pseudo-replication).
    xmap = {"control": 0, "PBC": 1}
    jitter = {"control": np.linspace(-0.06, 0.06, 5), "PBC": np.linspace(-0.06, 0.06, 5)}
    for condition in ["control", "PBC"]:
        g = composition[composition.condition == condition].sort_values("donor_id")
        color = BLUE if condition == "control" else VERMILLION
        xs = xmap[condition] + jitter[condition]
        ax_d.scatter(xs, 100 * g.lrNK_fraction, s=42, color=color, edgecolor="white", linewidth=0.65, zorder=3)
        mean = 100 * g.lrNK_fraction.mean()
        ax_d.plot([xmap[condition] - 0.18, xmap[condition] + 0.18], [mean, mean], color=INK, lw=2.0)
    p_comp = composition_test.loc[composition_test.metric == "lrNK_fraction", "p_value"].iloc[0]
    ax_d.set_xticks([0, 1])
    ax_d.set_xticklabels(["Control\n(n=5 donors)", "PBC\n(n=5 donors)"])
    ax_d.set_ylabel("lrNK / (lrNK + cNK), %")
    ax_d.set_xlim(-0.45, 1.45)
    ax_d.set_ylim(0, max(65, 100 * composition.lrNK_fraction.max() + 8))
    ax_d.set_title("NK composition: donor heterogeneity", loc="left", color=INK)
    ax_d.text(0.98, 0.95, f"Exact permutation P={p_comp:.3f}", transform=ax_d.transAxes, ha="right", va="top", fontsize=7.4, color=INK)
    ax_d.text(0.02, 0.03, "Each point is one donor", transform=ax_d.transAxes, fontsize=6.8, color=MUTED)
    style_axis(ax_d, "y")
    panel_label(ax_d, "D")

    # E: targeted donor-level pseudobulk effects, prioritising interpretable state changes.
    genes = ["IL12RB2", "TYK2", "STAT4", "IFNG"]
    cells = ["CD4_T", "cNK", "lrNK", "CD8_T"]
    ee = edge[(edge.scheme == "primary_conservative") & (edge.cell_filter == "all_final_singlets") &
              (edge.minimum_cells == 20) & edge.gene_symbol.isin(genes) & edge.cell_type.isin(cells)].copy()
    # TYK2 was not in the original target table; pull it from the matching all-gene files.
    if "TYK2" not in set(ee.gene_symbol):
        extra = []
        for cell in cells:
            path = HRA / "all_gene_results" / f"primary_conservative__all_final_singlets__{cell}_min20_edgeR.tsv"
            row = read(path)
            row = row[row.gene_symbol == "TYK2"].copy()
            row["cell_type"] = cell
            extra.append(row)
        ee = pd.concat([ee, *extra], ignore_index=True, sort=False)
    xpos = {g: i for i, g in enumerate(genes)}
    offsets = {"CD4_T": -0.24, "cNK": -0.08, "lrNK": 0.08, "CD8_T": 0.24}
    colors = {"CD4_T": GREEN, "cNK": ORANGE, "lrNK": VERMILLION, "CD8_T": BLUE}
    markers = {"CD4_T": "o", "cNK": "s", "lrNK": "D", "CD8_T": "^"}
    for cell in cells:
        sub = ee[ee.cell_type == cell]
        for _, row in sub.iterrows():
            low, high = ci_from_f(float(row.logFC), float(row.F))
            x = xpos[row.gene_symbol] + offsets[cell]
            ax_e.errorbar(x, row.logFC, yerr=[[row.logFC-low], [high-row.logFC]], fmt=markers[cell],
                          color=colors[cell], markersize=5.0, capsize=2.2, lw=1.1, zorder=3)
    ax_e.axhline(0, color="#6B7780", lw=0.9)
    ax_e.set_xticks(range(len(genes)))
    ax_e.set_xticklabels(genes, fontstyle="italic")
    ax_e.set_ylabel("PBC vs control log2 fold change")
    ax_e.set_title("Pseudobulk pathway effects", loc="left", color=INK)
    handles = [Line2D([0], [0], marker=markers[c], color="none", markerfacecolor=colors[c], markeredgecolor=colors[c], label=c.replace("_", " "), markersize=5.2) for c in cells]
    ax_e.legend(handles=handles, ncol=2, loc="upper right", columnspacing=0.9, handletextpad=0.3)
    ax_e.text(0.01, 0.02, "Exploratory: no target except CD8 T STAT4 passed full-transcriptome FDR 0.05", transform=ax_e.transAxes, fontsize=6.4, color=MUTED)
    style_axis(ax_e, "y")
    panel_label(ax_e, "E")

    # F: three complementary robustness metrics for IL12RB2.
    cells_f = ["CD4_T", "cNK", "lrNK", "CD8_T"]
    labels_f = ["CD4 T", "cNK", "lrNK", "CD8 T"]
    e_il = ee[ee.gene_symbol == "IL12RB2"].set_index("cell_type")
    p_il = permutation[(permutation.scheme == "primary_conservative") &
                       (permutation.cell_filter == "all_final_singlets") &
                       (permutation.minimum_cells == 20) &
                       (permutation.metric == "logCPM") &
                       (permutation.gene_symbol == "IL12RB2")].set_index("cell_type")
    l_il = lodo[(lodo.gene_symbol == "IL12RB2") & lodo.cell_type.isin(cells_f)].groupby("cell_type").first()
    y = np.arange(4)[::-1]
    for yi, cell in zip(y, cells_f):
        beta = float(e_il.loc[cell, "logFC"])
        low, high = ci_from_f(beta, float(e_il.loc[cell, "F"]))
        ax_f.errorbar(beta, yi, xerr=[[beta-low], [high-beta]], fmt="o", color=colors[cell], markersize=5.5, capsize=2.5)
        perm_p = float(p_il.loc[cell, "exact_PValue"])
        direction = float(l_il.loc[cell, "direction_preserved_fraction"])
        ax_f.text(3.55, yi, f"Pperm={perm_p:.3f}  LODO={direction*100:.0f}%", va="center", fontsize=6.3, color=INK)
    ax_f.axvline(0, color="#6B7780", lw=0.9)
    ax_f.set_yticks(y)
    ax_f.set_yticklabels(labels_f)
    ax_f.set_xlim(-1.2, 6.3)
    ax_f.set_xlabel("IL12RB2 robust QL log2 fold change")
    ax_f.set_title("IL12RB2 robustness analyses", loc="left", color=INK)
    ax_f.text(0.99, 0.02, "Direction stability ≠ external replication", transform=ax_f.transAxes, ha="right", fontsize=6.3, color=MUTED)
    style_axis(ax_f, "x")
    panel_label(ax_f, "F")

    fig.suptitle("Figure 5 | Liver single-cell localization of the IL12RB2–TYK2–STAT4–IFNG axis",
                 x=0.045, y=0.993, ha="left", fontsize=15.6, fontweight="bold", color=INK)
    fig.text(0.995, 0.005,
             "HRA008003: SingleR annotation and visualization at cell level; all disease contrasts and composition tests use donors as independent units.",
             ha="right", fontsize=7.1, color=MUTED)
    fig.savefig(OUT / "Fig05_single_cell_localization_final.pdf")
    fig.savefig(OUT / "Fig05_single_cell_localization_final.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
