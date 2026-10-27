#!/usr/bin/env python3
"""Generate review version of Supplementary Figure S1.

Screen-wide signal-level colocalization, prior sensitivity and SuSiE QC.
All plotted values are read from repository-relative locked analysis tables.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "data" / "derived" / "analysis" / "genetics" / "discovery"
OUT = REPO / "results" / "figures" / "supplementary"

P_COMBO = RESULTS / "23a_coloc_susie_combination_summary.tsv"
P_LOCUS = RESULTS / "23b_pbc_locus_finemapping_summary.tsv"
P_L10 = RESULTS / "20_onek1k_susie_concordance.tsv"
P_L3 = RESULTS / "20_onek1k_susie_concordance_l3_review.tsv"
for p in (P_COMBO, P_LOCUS, P_L10, P_L3):
    if not p.exists():
        raise FileNotFoundError(p)

combo = pd.read_csv(P_COMBO, sep="\t")
locus = pd.read_csv(P_LOCUS, sep="\t")
l10 = pd.read_csv(P_L10, sep="\t")
l3 = pd.read_csv(P_L3, sep="\t")

assert len(combo) == 23, f"Expected 23 combinations, found {len(combo)}"
assert len(locus) == 10, f"Expected 10 loci, found {len(locus)}"
assert locus["converged"].astype(str).str.lower().eq("true").all()

NAVY = "#17324D"
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
AMBER = "#E69F00"
VERMILION = "#D55E00"
PURPLE = "#8C6BB1"
RED = "#B94A48"
INK = "#26343F"
MID = "#687681"
LIGHT = "#E4E9ED"
PALE = "#F5F8F9"
WHITE = "#FFFFFF"

CLASS_COLOR = {
    "strong": GREEN,
    "candidate": BLUE,
    "exploratory": AMBER,
    "no_support": "#9AA4AC",
    "excluded_nonconvergent_eqtl": RED,
}
CLASS_LABEL = {
    "strong": "Strong",
    "candidate": "Candidate",
    "exploratory": "Exploratory",
    "no_support": "No support",
    "excluded_nonconvergent_eqtl": "Excluded: non-convergent",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8.0,
    "axes.titlesize": 10.0,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.3,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "legend.fontsize": 7.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def panel_header(ax, letter, title, subtitle=""):
    ax.text(-0.015, 1.055, letter, transform=ax.transAxes, ha="left", va="top",
            fontsize=14, fontweight="bold", color=NAVY)
    ax.text(0.043, 1.05, title, transform=ax.transAxes, ha="left", va="top",
            fontsize=10.1, fontweight="bold", color=INK)
    if subtitle:
        ax.text(0.043, 1.005, subtitle, transform=ax.transAxes, ha="left", va="top",
                fontsize=7.0, color=MID)


def rounded_box(ax, x, y, w, h, fc=WHITE, ec=LIGHT, lw=1.0):
    patch = FancyBboxPatch((x, y), w, h, transform=ax.transAxes,
                           boxstyle="round,pad=0.009,rounding_size=0.018",
                           facecolor=fc, edgecolor=ec, linewidth=lw, clip_on=False)
    ax.add_patch(patch)
    return patch


def bool_series(s):
    return s.astype(str).str.lower().eq("true")


# Display order: strict evidence class, then strict PP.H4; excluded remains last.
class_rank = {"strong": 0, "candidate": 1, "exploratory": 2,
              "no_support": 3, "excluded_nonconvergent_eqtl": 4}
combo["class_rank"] = combo["final_evidence_class"].map(class_rank)
combo = combo.sort_values(["class_rank", "PP.H4"], ascending=[True, False]).reset_index(drop=True)
combo["display"] = combo["gene_symbol"] + "  |  " + combo["sample_group"].str.replace("_", " ")

fig = plt.figure(figsize=(15.6, 11.4), facecolor=WHITE)
outer = fig.add_gridspec(2, 1, height_ratios=[1.57, 0.82], hspace=0.25,
                         left=0.075, right=0.978, top=0.91, bottom=0.065)
top = outer[0].subgridspec(1, 2, width_ratios=[1.22, 0.78], wspace=0.24)
bottom = outer[1].subgridspec(1, 2, width_ratios=[1.18, 0.82], wspace=0.26)
axA = fig.add_subplot(top[0, 0])
axB = fig.add_subplot(top[0, 1])
axC = fig.add_subplot(bottom[0, 0])
axD = fig.add_subplot(bottom[0, 1])

fig.suptitle("Screen-wide signal-level colocalization and model quality control",
             fontsize=15.2, fontweight="bold", color=NAVY, y=0.978)
fig.text(0.5, 0.951,
         "All 23 prespecified gene–cell combinations are shown; no result was removed because of statistical direction",
         ha="center", va="top", fontsize=8.5, color=MID)


# A | Complete strict-prior screen
y = np.arange(len(combo))
for yy, (_, r) in zip(y, combo.iterrows()):
    cls = r["final_evidence_class"]
    color = CLASS_COLOR[cls]
    if cls == "excluded_nonconvergent_eqtl":
        axA.hlines(yy, 0, r["PP.H4"], color=LIGHT, lw=1.1)
        axA.scatter(r["PP.H4"], yy, s=58, marker="x", color=color, lw=1.8, zorder=3)
        axA.text(1.015, yy, "not inferential", va="center", ha="left", fontsize=6.4,
                 color=color, fontweight="bold")
    else:
        axA.hlines(yy, 0, r["PP.H4"], color=color, alpha=0.42, lw=1.5)
        size = 34 + 3.0 * np.sqrt(max(float(r["credible_set_overlap_n"]), 0))
        axA.scatter(r["PP.H4"], yy, s=size, color=color, edgecolor=WHITE, lw=0.65, zorder=3)
        axA.text(1.015, yy, f"{r['PP.H4']:.3f}  |  overlap {int(r['credible_set_overlap_n'])}",
                 va="center", ha="left", fontsize=6.2, color=INK)

axA.axvspan(0.0, 0.1, color="#EEF1F3", zorder=0)
axA.axvspan(0.1, 0.5, color="#FFF7E5", zorder=0)
axA.axvspan(0.5, 0.9, color="#EBF4FA", zorder=0)
axA.axvspan(0.9, 1.0, color="#E7F5F0", zorder=0)
for xx in (0.1, 0.5, 0.9):
    axA.axvline(xx, color=WHITE, lw=1.3, zorder=1)
axA.set_xlim(0, 1.23)
axA.set_ylim(-0.8, len(combo)-0.2)
axA.set_yticks(y)
axA.set_yticklabels(combo["display"])
axA.invert_yaxis()
axA.set_xlabel(r"Signal-pair PP.H4 at strict prior $p_{12}=1\times10^{-6}$")
axA.grid(axis="x", color=LIGHT, lw=0.6, alpha=0.65)
axA.spines[["top", "right", "left"]].set_visible(False)
axA.tick_params(axis="y", length=0)
panel_header(axA, "A", "Complete prespecified discovery screen",
             "Marker size reflects credible-set overlap; the NFKB1 estimate is displayed only for audit")

# class separators and legend
for i in range(1, len(combo)):
    if combo.loc[i, "class_rank"] != combo.loc[i-1, "class_rank"]:
        axA.axhline(i-0.5, color="#C8D1D7", lw=0.8)
handles = []
for cls in ["strong", "candidate", "exploratory", "no_support", "excluded_nonconvergent_eqtl"]:
    marker = "x" if cls.startswith("excluded") else "o"
    handles.append(plt.Line2D([0], [0], marker=marker, ls="", color=CLASS_COLOR[cls],
                              markerfacecolor=CLASS_COLOR[cls], markersize=5.8,
                              label=CLASS_LABEL[cls]))
axA.legend(handles=handles, ncol=5, loc="lower left", bbox_to_anchor=(0.0, -0.13),
           frameon=False, columnspacing=1.2, handletextpad=0.35)


# B | Prior sensitivity for the same ordered screen
prior_cols = ["best_pph4_p12_1e6", "best_pph4_p12_5e6", "best_pph4_p12_1e5"]
heat = combo[prior_cols].to_numpy(float)
invalid_rows = combo["final_evidence_class"].eq("excluded_nonconvergent_eqtl").to_numpy()
heat_display = heat.copy()
heat_display[invalid_rows, :] = np.nan
cmap = LinearSegmentedColormap.from_list("pph4", ["#EEF1F3", "#F2C66D", "#56B4E9", "#009E73"])
cmap.set_bad("#E5E9EC")
im = axB.imshow(heat_display, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
axB.set_xticks([0, 1, 2])
axB.set_xticklabels([r"$1\times10^{-6}$"+"\nstrict", r"$5\times10^{-6}$",
                     r"$1\times10^{-5}$"+"\nrelaxed"])
axB.set_yticks(y)
axB.set_yticklabels(combo["analysis_id"] + "  " + combo["gene_symbol"], fontsize=6.4)
axB.tick_params(length=0)
for i in range(len(combo)+1):
    axB.axhline(i-0.5, color=WHITE, lw=0.9)
for j in range(4):
    axB.axvline(j-0.5, color=WHITE, lw=1.2)
for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        excluded = combo.loc[i, "final_evidence_class"] == "excluded_nonconvergent_eqtl"
        label = "X" if excluded else f"{heat[i, j]:.2f}"
        color = RED if excluded else (WHITE if heat[i, j] >= 0.68 else INK)
        axB.text(j, i, label, ha="center", va="center", fontsize=6.0,
                 color=color, fontweight="bold" if excluded else "normal")
for sp in axB.spines.values():
    sp.set_visible(False)
panel_header(axB, "B", "Prior sensitivity across all combinations",
             "Relaxed priors can upgrade candidates; final classes use the strict prior")
cbar = fig.colorbar(im, ax=axB, orientation="horizontal", fraction=0.035, pad=0.16, aspect=35)
cbar.set_label("PP.H4", fontsize=7.3)
cbar.ax.tick_params(labelsize=6.5, length=2)
axB.text(0.0, -0.315,
         "NFKB1 values arise from a non-convergent eQTL model and are crossed out rather than interpreted.",
         transform=axB.transAxes, fontsize=6.8, color=RED, va="top")


# C | Disease-side PBC SuSiE QC
locus = locus.sort_values("top_pip", ascending=True).reset_index(drop=True)
yc = np.arange(len(locus))
purity_norm = Normalize(vmin=0.5, vmax=1.0)
purity_cmap = LinearSegmentedColormap.from_list("purity", ["#E69F00", "#56B4E9", "#009E73"])
sizes = 55 + 48 * locus["credible_set_n"].astype(float)
sc = axC.scatter(locus["top_pip"], yc, s=sizes,
                 c=locus["minimum_cs_purity"], cmap=purity_cmap, norm=purity_norm,
                 edgecolor=WHITE, lw=0.8, zorder=3)
for yy, (_, r) in zip(yc, locus.iterrows()):
    axC.hlines(yy, 0, r["top_pip"], color=LIGHT, lw=1.2, zorder=1)
    axC.text(1.015, yy,
             f"CS={int(r['credible_set_n'])}; size {str(r['credible_set_sizes']).replace(';', '/')}; "
             f"purity={r['minimum_cs_purity']:.2f}",
             va="center", ha="left", fontsize=6.4, color=INK)
axC.set_yticks(yc)
axC.set_yticklabels(locus["pbc_locus"] + "  |  " + locus["lead_rsid"])
axC.set_xlim(0, 1.44)
axC.set_xlabel("Maximum disease SuSiE PIP (L=10)")
axC.grid(axis="x", color=LIGHT, lw=0.65)
axC.spines[["top", "right", "left"]].set_visible(False)
axC.tick_params(axis="y", length=0)
panel_header(axC, "C", "Disease-side fine-mapping QC across 10 loci",
             "All disease models converged; marker size denotes the number of 95% credible sets")
cbar2 = fig.colorbar(sc, ax=axC, orientation="horizontal", fraction=0.05, pad=0.17, aspect=35)
cbar2.set_label("Minimum credible-set purity", fontsize=7.2)
cbar2.ax.tick_params(labelsize=6.5, length=2)


# D | Convergence and exception handling
axD.set_axis_off()
panel_header(axD, "D", "Convergence audit and prespecified exception handling",
             "A failed model is excluded unless a lower-complexity review converges")

pbc_l10_ok = int(bool_series(l10["pbc_converged"]).sum())
eqtl_l10_ok = int(bool_series(l10["eqtl_converged"]).sum())
l3_ok = dict(zip(l3["analysis_id"], bool_series(l3["eqtl_converged"])))
final_ok = int((combo["final_evidence_class"] != "excluded_nonconvergent_eqtl").sum())

cards = [
    ("Disease SuSiE", "20 / 20", "L=3 and L=10 models", GREEN),
    ("eQTL SuSiE L=10", f"{eqtl_l10_ok} / 23", "two models reviewed at L=3", AMBER),
    ("L=3 review", f"{sum(l3_ok.values())} / {len(l3_ok)}", "IL12RB2 rescued; NFKB1 failed", BLUE),
    ("Final usable", f"{final_ok} / 23", "one model excluded", GREEN),
]
for i, (title, number, note, color) in enumerate(cards):
    x = 0.02 + (i % 2) * 0.49
    yy = 0.67 - (i // 2) * 0.25
    rounded_box(axD, x, yy, 0.46, 0.18, fc=PALE, ec=LIGHT)
    axD.text(x+0.025, yy+0.132, title, transform=axD.transAxes, fontsize=7.4,
             color=INK, fontweight="bold", va="center")
    axD.text(x+0.435, yy+0.132, number, transform=axD.transAxes, fontsize=10.2,
             color=color, fontweight="bold", ha="right", va="center")
    axD.text(x+0.025, yy+0.066, note, transform=axD.transAxes, fontsize=6.4,
             color=MID, va="center")

rounded_box(axD, 0.02, 0.08, 0.96, 0.27, fc=WHITE, ec=LIGHT)
axD.text(0.045, 0.298, "SA09  IL12RB2 / NK", transform=axD.transAxes,
         fontsize=7.5, color=INK, fontweight="bold")
axD.text(0.39, 0.298, "L10 failed  →  L3 converged  →  INCLUDED",
         transform=axD.transAxes, fontsize=7.2, color=GREEN, fontweight="bold")
axD.text(0.045, 0.225, "SA19  NFKB1 / CD4 TCM", transform=axD.transAxes,
         fontsize=7.5, color=INK, fontweight="bold")
axD.text(0.39, 0.225, "L10 failed  →  L3 failed  →  EXCLUDED",
         transform=axD.transAxes, fontsize=7.2, color=RED, fontweight="bold")
axD.text(0.045, 0.14,
         "Interpretation lock: PP.H4 from a non-convergent model is an audit trace, not biological evidence.",
         transform=axD.transAxes, fontsize=6.8, color=RED, fontweight="bold")


fig.text(0.075, 0.022,
         r"Evidence classes at strict $p_{12}=1\times10^{-6}$: strong PP.H4$\geq$0.90; candidate 0.50–0.90; exploratory 0.10–0.50; "
         "no support <0.10. Combinations sharing one PBC locus are not independent discoveries.",
         ha="left", va="bottom", fontsize=6.6, color=MID)

png = OUT / "FigS01_genetic_screen_coloc_QC_final.png"
pdf = OUT / "FigS01_genetic_screen_coloc_QC_final.pdf"
fig.savefig(pdf, bbox_inches="tight", facecolor=WHITE)
fig.savefig(png, dpi=300, bbox_inches="tight", facecolor=WHITE)
print(png)
print(pdf)
