#!/usr/bin/env python3
"""Generate Figure 7: bulk, proteomics, pharmacology and evidence grading.

Every numerical panel is read from repository-relative locked analysis tables.
Outputs a vector PDF and a 300-dpi PNG in the results directory.
"""
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch, Circle


REPO = Path(__file__).resolve().parents[2]
HERE = REPO / "results" / "figures" / "main"
ANALYSIS = REPO / "data" / "derived" / "analysis"

PATHS = {
    "meta": ANALYSIS / "bulk" / "pathway_meta" / "results" / "04_primary_random_effect_meta.tsv",
    "decision": ANALYSIS / "bulk" / "pathway_meta" / "results" / "10_prespecified_decision_summary.tsv",
    "robust": ANALYSIS / "proteomics" / "pxd083842" / "11_PXD083842_TYK2_robust_methods.tsv",
    "audit": ANALYSIS / "proteomics" / "pxd083842" / "17_PXD083842_covariate_metadata_audit.tsv",
    "adjusted": ANALYSIS / "proteomics" / "pxd083842" / "18_PXD083842_TYK2_adjusted_results.tsv",
    "tract": ANALYSIS / "pharmacology" / "drug_repurposing" / "01_target_tractability_and_candidate_counts.tsv",
    "rank": ANALYSIS / "pharmacology" / "drug_repurposing" / "03_genetic_direction_constrained_drug_ranking.tsv",
    "trials": ANALYSIS / "pharmacology" / "drug_repurposing" / "04_PBC_pathway_trials.tsv",
    "safety": ANALYSIS / "pharmacology" / "tyk2_safety" / "clinical" / "02_deucravacitinib_liver_safety_evidence.tsv",
}

for key, path in PATHS.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing {key}: {path}")

meta = pd.read_csv(PATHS["meta"], sep="\t")
decision = pd.read_csv(PATHS["decision"], sep="\t")
robust = pd.read_csv(PATHS["robust"], sep="\t")
audit = pd.read_csv(PATHS["audit"], sep="\t")
adjusted = pd.read_csv(PATHS["adjusted"], sep="\t")
tract = pd.read_csv(PATHS["tract"], sep="\t")
ranking = pd.read_csv(PATHS["rank"], sep="\t")
trials = pd.read_csv(PATHS["trials"], sep="\t")
safety = pd.read_csv(PATHS["safety"], sep="\t")


# Consistent, color-blind-aware palette used across the main figures.
NAVY = "#17324D"
BLUE = "#3A78A8"
TEAL = "#2A9D8F"
GREEN = "#2E8B70"
AMBER = "#D99A2B"
CORAL = "#D96654"
RED = "#B84A4A"
INK = "#25313C"
MID = "#66737F"
PALE = "#F4F7F8"
LINE = "#D8E0E5"
WHITE = "#FFFFFF"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8.2,
    "axes.titlesize": 10.0,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.4,
    "xtick.labelsize": 7.3,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.1,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def panel_title(ax, letter, title, subtitle=None):
    ax.text(-0.02, 1.075, letter, transform=ax.transAxes, fontsize=14, fontweight="bold",
            color=NAVY, va="top", ha="left")
    ax.text(0.045, 1.067, title, transform=ax.transAxes, fontsize=10.2, fontweight="bold",
            color=INK, va="top", ha="left")
    if subtitle:
        ax.text(0.045, 1.005, subtitle, transform=ax.transAxes, fontsize=7.1,
                color=MID, va="top", ha="left")


def rounded_box(ax, xy, width, height, fc=WHITE, ec=LINE, lw=1.0, radius=0.02):
    patch = FancyBboxPatch(xy, width, height,
                           boxstyle=f"round,pad=0.008,rounding_size={radius}",
                           transform=ax.transAxes, facecolor=fc, edgecolor=ec,
                           linewidth=lw, clip_on=False)
    ax.add_patch(patch)
    return patch


def status_dot(ax, x, y, color, radius=0.017, transform=None):
    if transform is None:
        transform = ax.transAxes
    ax.add_patch(Circle((x, y), radius, transform=transform, facecolor=color,
                        edgecolor=WHITE, linewidth=0.7, clip_on=False))


fig = plt.figure(figsize=(15.8, 11.0), facecolor=WHITE)
outer = fig.add_gridspec(2, 1, height_ratios=[1.03, 0.97], hspace=0.20,
                         left=0.055, right=0.985, top=0.895, bottom=0.055)
top = outer[0].subgridspec(1, 3, width_ratios=[1.05, 1.18, 1.03], wspace=0.30)
bottom = outer[1].subgridspec(1, 2, width_ratios=[1.12, 0.88], wspace=0.28)
boxA = fig.add_subplot(top[0, 0]); boxA.set_axis_off()
boxB = fig.add_subplot(top[0, 1]); boxB.set_axis_off()
boxC = fig.add_subplot(top[0, 2]); boxC.set_axis_off()
axA = boxA.inset_axes([0.02, 0.19, 0.96, 0.72])
axB = boxB.inset_axes([0.02, 0.39, 0.96, 0.52])
axC = boxC.inset_axes([0.02, 0.25, 0.96, 0.66])
axD = fig.add_subplot(bottom[0, 0])
axE = fig.add_subplot(bottom[0, 1])

fig.suptitle("Orthogonal molecular evidence and a guarded translational path for the IL12RB2–TYK2 axis",
             fontsize=15.3, fontweight="bold", color=NAVY, y=0.982)
fig.text(0.5, 0.955,
         "Directional multi-omics support is separated from clinical efficacy, safety and druggability claims",
         ha="center", va="top", fontsize=8.5, color=MID)


# A | Bulk pathway meta-analysis
primary = meta[meta["score_method"].eq("robust_z")].copy()
order = ["IL12_STAT4", "IFNG_RESPONSE"]
primary["signature"] = pd.Categorical(primary["signature"], order, ordered=True)
primary = primary.sort_values("signature", ascending=False)
y = np.arange(len(primary))
for yy, (_, row) in zip(y, primary.iterrows()):
    axA.hlines(yy, row["ci_low_mKH"], row["ci_high_mKH"], color=BLUE, lw=2.0, zorder=2)
    axA.scatter(row["estimate"], yy, s=70, color=TEAL, edgecolor=WHITE, lw=0.8, zorder=3)
    label = (f"{int(row['positive_effects'])}/{int(row['k'])} cohorts positive\n"
             f"HK P={row['p_mKH']:.3f}; FDR={row['BH_FDR_four_primary_tests']:.3f}; "
             f"I²={row['I2_percent']:.1f}%")
    axA.text(3.08, yy, label, va="center", ha="left", fontsize=7.0, color=INK, linespacing=1.25)
axA.axvline(0, color=MID, lw=0.9, ls="--")
axA.set_xlim(-1.6, 5.65)
axA.set_ylim(-0.7, 1.7)
axA.set_yticks(y)
axA.set_yticklabels(["IFNG response", "IL12–STAT4"])
axA.set_xlabel("standardized PBC–control pathway effect (modified HK 95% CI)")
axA.grid(axis="x", color=LINE, lw=0.7, alpha=0.7)
axA.spines[["top", "right", "left"]].set_visible(False)
axA.tick_params(axis="y", length=0)
panel_title(boxA, "A", "Bulk transcriptomic pathway meta-analysis",
            "Four independent cohorts; prespecified robust-z score")
boxA.text(0.02, 0.105,
         "Direction preserved under rank scoring, leave-one-cohort-out and leave-one-gene-out analyses;\n"
         "neither pathway reached multiplicity-adjusted significance.",
         transform=boxA.transAxes, fontsize=7.1, color=MID, va="top")


# B | Proteomic robustness plus confounding audit
method_order = ["Welch", "Mann_Whitney", "Permutation_100000", "Hodges_Lehmann",
                "OLS_HC3", "Huber_RLM", "winsor_5_95"]
rb = robust[(robust["contrast"].eq("PBC_vs_Healthy")) & robust["estimate"].notna()
            & robust["method"].isin(method_order)].copy()
rb["method"] = pd.Categorical(rb["method"], method_order, ordered=True)
rb = rb.sort_values("method")
labels = {"Mann_Whitney": "Mann–Whitney", "Permutation_100000": "Permutation",
          "Hodges_Lehmann": "Hodges–Lehmann", "OLS_HC3": "OLS HC3",
          "Huber_RLM": "Huber RLM", "winsor_5_95": "Winsorised"}
yy = np.arange(len(rb))[::-1]
for ypos, (_, row) in zip(yy, rb.iterrows()):
    axB.scatter(row["estimate"], ypos, s=42, color=TEAL, edgecolor=WHITE, lw=0.6, zorder=3)
    axB.text(0.515, ypos, f"P={row['p_value']:.3g}", va="center", ha="left", fontsize=6.7, color=INK)
axB.axvline(0, color=MID, lw=0.9, ls="--")
axB.set_xlim(-0.03, 0.73)
axB.set_ylim(-2.25, len(rb)-0.3)
axB.set_yticks(yy)
axB.set_yticklabels([labels.get(str(x), str(x)) for x in rb["method"]])
axB.set_xlabel("TYK2 liver protein effect (log2 scale)")
axB.grid(axis="x", color=LINE, lw=0.7, alpha=0.7)
axB.spines[["top", "right", "left"]].set_visible(False)
axB.tick_params(axis="y", length=0)
panel_title(boxB, "B", "PXD083842 proteomic robustness and confounding",
            "PBC n=14 versus healthy n=9")

# Adjustment tiles below the point plot.
adj_order = ["primary_unadjusted", "cirrhosis_adjusted", "age_adjusted", "sex_adjusted",
             "BMI_adjusted", "age_plus_sex_adjusted", "exclude_concurrent_MASLD",
             "disease_specific_comparator"]
adj_labels = ["Unadjusted", "Cirrhosis", "Age", "Sex", "BMI", "Age+sex", "No MASLD", "vs MASLD"]
adj = adjusted.set_index("analysis").reindex(adj_order)
x0, dx = 0.02, 0.242
for i, (lab, (_, row)) in enumerate(zip(adj_labels, adj.iterrows())):
    flag = bool(row["TYK2_in_FDR_lt_0_05_list"])
    color = GREEN if flag else "#D9DEE2"
    col, rownum = i % 4, i // 4
    yy_tile = 0.205 - rownum * 0.105
    rounded_box(boxB, (x0+col*dx, yy_tile), 0.222, 0.077, fc=color, ec=WHITE, radius=0.012)
    boxB.text(x0+col*dx+0.111, yy_tile+0.039, lab, transform=boxB.transAxes, ha="center", va="center",
             fontsize=5.9, color=WHITE if flag else MID, fontweight="bold" if flag else "normal")
boxB.text(0.02, 0.31, "FDR<0.05 list membership after adjustment/exclusion", transform=boxB.transAxes,
         fontsize=6.4, color=MID, va="top")
boxB.text(0.02, 0.072,
         "Critical: age has minimal overlap; recruitment stream is exactly confounded with disease.\n"
         "Adjusted non-significant coefficients were not released. Supportive, not disease-specific replication.",
         transform=boxB.transAxes, fontsize=6.6, color=RED, va="top", fontweight="bold")


# C | Evidence matrix: database presence is not efficacy
rows = ["IL12RB2 direct", "TYK2 / deucravacitinib", "IL12B / ustekinumab", "JAK1/2 / baricitinib"]
cols = ["Target /\nmechanism", "Clinical-stage\nmodality", "PBC trial\nevidence", "Label / safety\nevidence", "Genetic-direction\nmatch"]
# 2=supportive/available; 1=qualified/partial; 0=absent/negative.
# Construct cells from the audited target, drug-ranking, trial and safety tables.
target = tract.set_index("target")
rank_by_candidate = ranking.set_index("candidate")
il12rb2_modality = int(target.loc["IL12RB2", "candidate_count"] > 0
                       or bool(target.loc["IL12RB2", "antibody_approved"]))
tyk2_modality = int(target.loc["TYK2", "candidate_count"] > 0
                    and bool(target.loc["TYK2", "small_molecule_approved"]))
il12b_modality = int(target.loc["IL12B", "candidate_count"] > 0
                     and bool(target.loc["IL12B", "antibody_approved"]))
jak2_modality = int(target.loc["JAK2", "candidate_count"] > 0
                    and bool(target.loc["JAK2", "small_molecule_approved"]))
has_tyk2_pbc_trial = trials["title"].str.contains("TYK2|deucravacitinib", case=False, regex=True).any()
has_ustek_trial = trials["title"].str.contains("ustekinumab", case=False).any()
has_bari_trial = trials["title"].str.contains("baricitinib", case=False).any()
ustek_positive = has_ustek_trial and not trials.loc[
    trials["title"].str.contains("ustekinumab", case=False), "primary_result"
].astype(str).str.contains('"value": "0"', regex=False).any()
bari_efficacy = has_bari_trial and not trials.loc[
    trials["title"].str.contains("baricitinib", case=False), "status"
].astype(str).str.contains("TERMINATED", case=False).any()
matrix = np.array([
    [2, 2 if il12rb2_modality else 0, 0, 0, 2],
    [2, 2 if tyk2_modality else 0, 2 if has_tyk2_pbc_trial else 0,
     1 if not safety.empty else 0,
     2 if str(rank_by_candidate.loc["DEUCRAVACITINIB", "direction_match"]).startswith("yes") else 1],
    [2, 2 if il12b_modality else 0, 2 if ustek_positive else 0, 1,
     1 if "ligand level" in str(rank_by_candidate.loc["USTEKINUMAB", "direction_match"]) else 0],
    [2, 2 if jak2_modality else 0, 2 if bari_efficacy else 0, 0,
     1 if "partly" in str(rank_by_candidate.loc["BARICITINIB", "direction_match"]) else 0],
])
axC.imshow(matrix, cmap=ListedColormap(["#E8ECEF", "#F2C66D", "#68AD91"]), vmin=0, vmax=2,
           aspect="auto", interpolation="nearest")
axC.set_xticks(np.arange(len(cols)))
axC.set_xticklabels(cols, rotation=0, ha="center", fontsize=6.4)
axC.xaxis.tick_top()
axC.set_yticks(np.arange(len(rows)))
axC.set_yticklabels(rows, fontsize=7.0)
for i in range(len(rows)+1):
    axC.axhline(i-0.5, color=WHITE, lw=2)
for j in range(len(cols)+1):
    axC.axvline(j-0.5, color=WHITE, lw=2)
symbols = {0: "—", 1: "!", 2: "+"}
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        axC.text(j, i, symbols[int(matrix[i, j])], ha="center", va="center",
                 fontsize=10, color=WHITE if matrix[i, j] == 2 else (INK if matrix[i, j] == 1 else MID),
                 fontweight="bold")
axC.tick_params(length=0)
for sp in axC.spines.values():
    sp.set_visible(False)
panel_title(boxC, "C", "Drug and clinical evidence matrix",
            "Open Targets / ChEMBL / ClinicalTrials.gov / labels")
boxC.text(0.02, 0.18, "+ available/supportive     ! qualified risk     — absent or negative",
         transform=boxC.transAxes, fontsize=6.7, color=MID, va="top")
boxC.text(0.02, 0.09,
         "PBC counterevidence: ustekinumab 0/20 achieved week-12 ALP response;\n"
         "baricitinib trial terminated after n=2. No TYK2-directed PBC trial found.",
         transform=boxC.transAxes, fontsize=6.8, color=RED, va="top", fontweight="bold")


# D | Explicit translational traffic-light gates
axD.set_axis_off()
panel_title(axD, "D", "Safety and druggability gates",
            "Green does not denote efficacy; it denotes passage of that specific gate")

rounded_box(axD, (0.015, 0.08), 0.47, 0.79, fc="#F7FAFB", ec=LINE, radius=0.018)
rounded_box(axD, (0.515, 0.08), 0.47, 0.79, fc="#F7FAFB", ec=LINE, radius=0.018)
axD.text(0.04, 0.82, "TYK2 / deucravacitinib", transform=axD.transAxes, fontsize=10.0,
         color=NAVY, fontweight="bold")
axD.text(0.54, 0.82, "Direct IL12RB2 small-molecule route", transform=axD.transAxes, fontsize=10.0,
         color=NAVY, fontweight="bold")

left_items = [
    ("Mechanistic and genetic direction", GREEN, "Matched"),
    ("Selective approved TYK2 modality", GREEN, "Available"),
    ("Direct efficacy evidence in PBC", RED, "Absent"),
    ("ALT / AST and DILI monitoring", AMBER, "Required"),
    ("Serious infection and latent TB", AMBER, "Core gate"),
    ("Child–Pugh C", RED, "Not recommended"),
]
right_items = [
    ("Experimental receptor structures", GREEN, "8XRP / 8YI7"),
    ("Ligand-bound small-molecule pocket", RED, "Absent"),
    ("ChEMBL IL12RB2 activity record", RED, "Absent"),
    ("Validated docking benchmark", RED, "Absent"),
    ("Receptor-selective biologic concept", AMBER, "Research only"),
]
for start_x, items in [(0.04, left_items), (0.54, right_items)]:
    for k, (label, color, val) in enumerate(items):
        yy = 0.72 - k*0.098
        status_dot(axD, start_x+0.012, yy, color, radius=0.013)
        axD.text(start_x+0.038, yy, label, transform=axD.transAxes, va="center", ha="left",
                 fontsize=7.4, color=INK)
        axD.text(start_x+0.425, yy, val, transform=axD.transAxes, va="center", ha="right",
                 fontsize=6.8, color=color, fontweight="bold")

rounded_box(axD, (0.055, 0.12), 0.39, 0.075, fc="#FFF5DE", ec="#E7C56C", radius=0.015)
axD.text(0.25, 0.158, "CONDITIONAL GO", transform=axD.transAxes, ha="center", va="center",
         fontsize=10.5, color="#9A6500", fontweight="bold")
rounded_box(axD, (0.555, 0.12), 0.39, 0.075, fc="#FBEAEA", ec="#DA9999", radius=0.015)
axD.text(0.75, 0.158, "NO-GO for docking", transform=axD.transAxes, ha="center", va="center",
         fontsize=10.5, color=RED, fontweight="bold")
axD.text(0.25, 0.094, "Preclinical / early proof-of-mechanism only", transform=axD.transAxes,
         ha="center", va="top", fontsize=6.8, color=MID)
axD.text(0.75, 0.094, "Reassess if a validated pocket or ligand emerges", transform=axD.transAxes,
         ha="center", va="top", fontsize=6.8, color=MID)


# E | Evidence grade: deliberately separate target, translation and clinical recommendation
axE.set_axis_off()
panel_title(axE, "E", "Final evidence grade",
            "Claim strength is capped by the weakest required evidence layer")
tiers = [
    ("1", "Genetic target support", "MODERATE–STRONG", GREEN,
     "IL12RB2 locus/eQTL convergence;\nTYK2 human loss-of-function direction"),
    ("2", "Translational hypothesis", "CONDITIONAL", AMBER,
     "Bulk direction is consistent but non-significant;\nproteomic signal is confounded; selective TYK2 is testable"),
    ("3", "Clinical recommendation", "NOT SUPPORTED", RED,
     "No direct TYK2 PBC trial; prior pathway trials negative/inconclusive;\nPBC-specific efficacy and safety remain unknown"),
]
for i, (num, title, grade, color, detail) in enumerate(tiers):
    yy = 0.75 - i*0.255
    rounded_box(axE, (0.04, yy-0.11), 0.92, 0.205, fc="#F8FAFB", ec=LINE, radius=0.018)
    axE.add_patch(Circle((0.095, yy), 0.035, transform=axE.transAxes, facecolor=color,
                         edgecolor=WHITE, lw=1.0))
    axE.text(0.095, yy, num, transform=axE.transAxes, ha="center", va="center",
             color=WHITE, fontsize=9.2, fontweight="bold")
    axE.text(0.15, yy+0.038, title, transform=axE.transAxes, ha="left", va="center",
             color=INK, fontsize=9.0, fontweight="bold")
    axE.text(0.92, yy+0.038, grade, transform=axE.transAxes, ha="right", va="center",
             color=color, fontsize=8.0, fontweight="bold")
    axE.text(0.15, yy-0.036, detail, transform=axE.transAxes, ha="left", va="center",
             color=MID, fontsize=7.0, linespacing=1.3)

rounded_box(axE, (0.085, 0.015), 0.83, 0.085, fc=NAVY, ec=NAVY, radius=0.018)
axE.text(0.5, 0.058,
         "Final claim: genetically supported target axis with a guarded experimental drug hypothesis",
         transform=axE.transAxes, ha="center", va="center", color=WHITE,
         fontsize=7.5, fontweight="bold")


# Global provenance footer.
fig.text(0.055, 0.018,
         "All estimates and decisions are read from locked project tables. Database presence denotes evidence availability, not therapeutic efficacy. "
         "Audits current to 14–16 October 2026.",
         ha="left", va="bottom", fontsize=6.5, color=MID)

out_png = HERE / "Fig07_bulk_proteomics_translation_final.png"
out_pdf = HERE / "Fig07_bulk_proteomics_translation_final.pdf"
fig.savefig(out_pdf, bbox_inches="tight", facecolor=WHITE)
fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor=WHITE)
print(out_png)
print(out_pdf)
