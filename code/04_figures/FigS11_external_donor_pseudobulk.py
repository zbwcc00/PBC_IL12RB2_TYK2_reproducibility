from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
HERE = REPO / "results" / "figures" / "supplementary"
SOURCE = REPO / "data" / "derived" / "figures" / "supplementary" / "FigS11_Andrews_external_donor_pseudobulk.tsv"

INK = "#263746"
PBC = "#2878B5"
NORMAL = "#8B949E"
GRID = "#DCE4EA"
ORANGE = "#D55E00"


def panel(ax, data, modality, target, letter, note, min_cells=10):
    q = data[(data.modality == modality) & (data.target == target) & (data.n_cells >= min_cells)].copy()
    groups = [("Normal", NORMAL), ("PBC", PBC)]
    rng = np.random.default_rng(20261020)
    for xpos, (group, color) in enumerate(groups):
        z = q[q.group == group].sort_values("donor_id")
        jitter = rng.uniform(-0.07, 0.07, len(z))
        ax.scatter(np.repeat(xpos, len(z)) + jitter, z.logCPM, s=54, color=color,
                   edgecolor="white", linewidth=0.7, zorder=3)
        if len(z):
            mean = z.logCPM.mean()
            ax.plot([xpos - 0.22, xpos + 0.22], [mean, mean], color=color, lw=2.0, zorder=2)
        offset_patterns = [(-20, 28), (-8, 42), (4, 28), (16, 42), (-14, 54), (10, 54)]
        for idx, (dx, (_, row)) in enumerate(zip(jitter, z.iterrows())):
            x_offset, y_offset = offset_patterns[idx % len(offset_patterns)]
            vertical_alignment = "bottom"
            if row.logCPM >= q.logCPM.max() - 0.05:
                y_offset = -25 - (idx % 2) * 14
                vertical_alignment = "top"
            ax.annotate(row.donor_id, (xpos + dx, row.logCPM), xytext=(x_offset, y_offset),
                        textcoords="offset points", fontsize=6.5, rotation=34, color=INK,
                        ha="right" if x_offset < 0 else "left", va=vertical_alignment,
                        rotation_mode="anchor")
    ax.set_xticks([0, 1], ["Normal", "PBC"])
    ax.set_xlim(-0.48, 1.48)
    ax.margins(y=0.13)
    ax.set_ylabel("IL12RB2 donor pseudobulk\nlog2(CPM + 0.5)")
    ax.set_title(f"{modality} | {target}", loc="left", color=INK, fontweight="bold", pad=25)
    ax.text(0.98, 1.01, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=7.1, color=INK,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=GRID, alpha=0.94))
    ax.text(-0.12, 1.14, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", color=INK)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    data = pd.read_csv(SOURCE, sep="\t")
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8.5, "figure.dpi": 300, "pdf.fonttype": 42,
                         "ps.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.2), constrained_layout=True)
    panel(axes[0, 0], data, "scRNA", "lrNK", "A",
          "PBC n=2; normal n=2 after ≥10-cell filter\nΔ=4.604; exact one-sided P=0.1667")
    panel(axes[0, 1], data, "scRNA", "T_NK_union", "B",
          "PBC n=2; normal n=6\nΔ=5.626; exact one-sided P=0.0357")
    panel(axes[1, 0], data, "snRNA", "lrNK", "C",
          "PBC n=2; normal n=3\nΔ=0.004; exact one-sided P=0.5000")
    panel(axes[1, 1], data, "snRNA", "T_NK_union", "D",
          "PBC n=2; normal n=3\nΔ=−0.226; exact one-sided P=0.6000")
    fig.suptitle("Supplementary Figure S11 | External Andrews liver-atlas donor-level PBC–normal audit",
                 x=0.02, ha="left", fontsize=14.3, fontweight="bold", color=INK)
    fig.text(0.99, -0.035,
             "Points are donors. The two modalities are technical arms of one study; PBC001 is shared and is not counted twice.",
             ha="right", fontsize=7.2, color=ORANGE)
    fig.savefig(HERE / "FigS11_Andrews_external_donor_pseudobulk_final.png", dpi=300, bbox_inches="tight")
    fig.savefig(HERE / "FigS11_Andrews_external_donor_pseudobulk_final.pdf", bbox_inches="tight")
    fig.savefig(HERE / "FigS11_Andrews_external_donor_pseudobulk_final.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
