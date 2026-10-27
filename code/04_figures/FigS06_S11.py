from __future__ import annotations

import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
ANALYSIS = REPO / "data" / "derived" / "analysis"
OUT = REPO / "results" / "figures" / "supplementary"
SRC = REPO / "data" / "derived" / "figures" / "supplementary"
OUT.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)

HRA = ANALYSIS / "single_cell" / "hra008003"
QC = HRA / "qc"
SOUPX = HRA / "soupx"
ANNO = HRA / "annotation"
PB = HRA / "pseudobulk"
BULK = ANALYSIS / "bulk" / "validation" / "results"
PATHWAY = ANALYSIS / "bulk" / "pathway_meta" / "results"
PXD = ANALYSIS / "proteomics" / "pxd083842"
TYK2 = ANALYSIS / "pharmacology" / "tyk2_safety" / "locus"

BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
YELLOW = "#F0E442"
VERMILLION = "#D55E00"
PINK = "#CC79A7"
GRAY = "#9AA1A8"
LIGHTGRAY = "#E6E9ED"
INK = "#263238"
NAVY = "#17324D"
COLORS = [BLUE, ORANGE, GREEN, VERMILLION, SKY, PINK, "#7B61A8", "#6C757D", "#A6761D", "#1B9E77"]


plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.2,
    "axes.titlesize": 10.2,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.8,
    "legend.fontsize": 7.2,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def style(ax, xgrid=False, ygrid=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if xgrid:
        ax.grid(axis="x", color="#D8DDE3", lw=0.6, alpha=0.65)
    if ygrid:
        ax.grid(axis="y", color="#D8DDE3", lw=0.6, alpha=0.65)
    ax.set_axisbelow(True)


def label(ax, letter):
    ax.text(-0.10, 1.08, letter, transform=ax.transAxes, fontsize=12.5,
            fontweight="bold", va="top", ha="center", color="white",
            bbox=dict(boxstyle="round,pad=0.25", facecolor=NAVY, edgecolor="none"))


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    plt.close(fig)


def sample_order(df):
    controls = sorted([x for x in df["sample_id"].unique() if str(x).startswith("control")])
    pbc = sorted([x for x in df["sample_id"].unique() if str(x).startswith("pbc")])
    return controls + pbc


def short_sample(x):
    return str(x).replace("control_liver", "Control ").replace("pbc_liver", "PBC ")


def load_qc():
    frames = [read_tsv(p) for p in sorted((QC / "sample_summary").glob("*_summary.tsv"))]
    q = pd.concat(frames, ignore_index=True)
    soup, target = [], []
    for d in sorted(SOUPX.iterdir()):
        if not d.is_dir():
            continue
        s = read_tsv(d / "SoupX_summary.tsv")
        g = read_tsv(d / "candidate_gene_sensitivity.tsv")
        g = g[g["gene"] == "IL12RB2"].copy()
        soup.append(s)
        target.append(g.assign(sample_id=d.name))
    return q, pd.concat(soup, ignore_index=True), pd.concat(target, ignore_index=True)


def fig_s1():
    q, soup, target = load_qc()
    order = sample_order(q)
    q["sample_id"] = pd.Categorical(q["sample_id"], order, ordered=True)
    q = q.sort_values("sample_id")
    merged = q.merge(soup, on="sample_id").merge(target, on="sample_id")
    merged.to_csv(SRC / "FigS6_HRA008003_QC_source.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.5, 9.0), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.98, bottom=0.09, top=0.84,
                          hspace=0.42, wspace=0.30)
    axa, axb, axc, axd = [fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    x = np.arange(len(q))
    stages = [
        ("Cell-called", "union_cell_called", GRAY),
        ("Core QC", "passes_core_qc", SKY),
        ("Final singlets", "final_singlets", GREEN),
    ]
    for name, col, color in stages:
        frac = q[col] / q["union_cell_called"] * 100
        axa.plot(x, frac, marker="o", ms=4.5, lw=1.5, label=name, color=color)
    axa.axvline(4.5, color="#B8BEC5", ls="--", lw=0.9)
    axa.set_xticks(x, [short_sample(v) for v in q["sample_id"].astype(str)], rotation=45, ha="right")
    axa.set_ylabel("Retained relative to cell-called (%)")
    axa.set_ylim(55, 103)
    axa.set_title("Cell-calling and quality-control retention", loc="left")
    axa.legend(ncol=3, loc="lower left")
    style(axa, ygrid=True); label(axa, "A")

    bar_colors = [BLUE if c == "control" else VERMILLION for c in q["condition"]]
    axb.bar(x, q["final_singlets"], color=bar_colors, width=0.72)
    for xi, value in zip(x, q["final_singlets"]):
        axb.text(xi, value + q["final_singlets"].max() * 0.018, f"{int(value):,}",
                 ha="center", va="bottom", fontsize=6.4)
    axb.axvline(4.5, color="#B8BEC5", ls="--", lw=0.9)
    axb.set_xticks(x, [short_sample(v) for v in q["sample_id"].astype(str)], rotation=45, ha="right")
    axb.set_ylabel("Final singlets")
    axb.set_ylim(0, q["final_singlets"].max() * 1.105)
    axb.set_title(f"Final donor-level cell yield (total={int(q['final_singlets'].sum()):,})", loc="left")
    axb.legend(handles=[Line2D([0], [0], color=BLUE, lw=7, label="Control"), Line2D([0], [0], color=VERMILLION, lw=7, label="PBC")])
    style(axb, ygrid=True); label(axb, "B")

    metrics = ["final_singlet_median_UMI", "final_singlet_median_genes", "final_singlet_median_percent_mt"]
    names = ["Median UMI", "Median genes", "Median mitochondrial %"]
    z = q[metrics].apply(lambda c: (c - c.mean()) / c.std(ddof=0))
    im = axc.imshow(z.T, aspect="auto", cmap=LinearSegmentedColormap.from_list("qc", [BLUE, "#FFFFFF", VERMILLION]), vmin=-2, vmax=2)
    axc.set_xticks(x, [short_sample(v) for v in q["sample_id"].astype(str)], rotation=45, ha="right")
    axc.set_yticks(np.arange(3), names)
    for iy, col in enumerate(metrics):
        for ix, value in enumerate(q[col]):
            text = f"{value:,.0f}" if col != "final_singlet_median_percent_mt" else f"{value:.1f}"
            axc.text(ix, iy, text, ha="center", va="center", fontsize=6.3, color=INK)
    fig.colorbar(im, ax=axc, fraction=0.035, pad=0.02, label="Row-standardized value")
    axc.set_title("Per-sample final-singlet QC metrics", loc="left")
    label(axc, "C")

    merged = merged.set_index("sample_id").loc[order].reset_index()
    axd.scatter(merged["median_contamination_fraction"] * 100, merged["UMI_retention_fraction"] * 100,
                s=55, c=[BLUE if x == "control" else VERMILLION for x in merged["condition"]], edgecolor="white", lw=0.6)
    for _, r in merged.iterrows():
        xx = r["median_contamination_fraction"] * 100
        yy = r["UMI_retention_fraction"] * 100
        dx, dy, ha = (3, 2, "left")
        if yy >= 99.85:
            dy = -10
        if xx >= 3.35:
            dx, ha = -4, "right"
        axd.annotate(short_sample(r["sample_id"]), (xx, yy), xytext=(dx, dy),
                     textcoords="offset points", fontsize=6.2, ha=ha)
    axd.axhline(100, color="#B8BEC5", lw=0.8)
    axd.set_xlabel("SoupX estimated contamination (%)")
    axd.set_ylabel("IL12RB2 UMI retained after correction (%)")
    axd.set_title("Ambient-RNA correction preserves target counts", loc="left")
    axd.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE,
                               markeredgecolor="white", markersize=6, label="Control"),
                         Line2D([0], [0], marker="o", color="none", markerfacecolor=VERMILLION,
                               markeredgecolor="white", markersize=6, label="PBC")],
                loc="lower left", bbox_to_anchor=(0.0, 0.12), handletextpad=0.3)
    axd.text(0.02, 0.05, f"Pooled retention={target['corrected_final_singlet_UMI'].sum()/target['raw_final_singlet_UMI'].sum()*100:.1f}%",
             transform=axd.transAxes, color=GREEN, fontweight="bold")
    style(axd, xgrid=True, ygrid=True); label(axd, "D")
    fig.suptitle("Supplementary Figure S6 | HRA008003 BAM reconstruction and quality control",
                 fontsize=14, fontweight="bold", x=0.08, y=0.975, ha="left", color=NAVY)
    fig.text(0.08, 0.938,
             "Donor-resolved reconstruction, cell retention and target-specific ambient-RNA sensitivity.",
             fontsize=8.5, color="#687681")
    save(fig, "FigS06_HRA008003_QC_final")


def fig_s2():
    counts = read_tsv(ANNO / "annotation_cell_counts.tsv")
    meta = read_tsv(ANNO / "all_processed_cell_metadata.tsv.gz", compression="gzip")
    order = sample_order(counts)
    label_totals = counts.groupby("conservative_label")["n_cells"].sum().sort_values(ascending=False)
    top = list(label_totals.head(9).index)
    counts["display_label"] = np.where(counts["conservative_label"].isin(top), counts["conservative_label"], "Other")
    comp = counts.groupby(["sample_id", "condition", "display_label"], as_index=False)["n_cells"].sum()
    piv = comp.pivot(index="sample_id", columns="display_label", values="n_cells").fillna(0).reindex(order)
    prop = piv.div(piv.sum(axis=1), axis=0) * 100
    ann = meta.groupby(["conservative_label", "condition"], as_index=False).agg(
        n_cells=("cell_id", "size"), median_delta=("annotation_delta", "median"),
        pruned_fraction=("annotation_pruned", "mean"), median_umi=("total_UMI_SoupX", "median"))
    tnk = meta[meta["conservative_label"].isin(["lrNK", "cNK", "CD4_T", "CD8_T", "TNK_mixed", "T_other"])].copy()
    cross = pd.crosstab(tnk["conservative_label"], tnk["tnk_refined_label"], normalize="index") * 100
    comp.to_csv(SRC / "FigS7_annotation_composition.tsv", sep="\t", index=False)
    ann.to_csv(SRC / "FigS7_annotation_confidence.tsv", sep="\t", index=False)
    cross.to_csv(SRC / "FigS7_TNK_refinement_concordance.tsv", sep="\t")

    fig = plt.figure(figsize=(13.5, 9.8), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.98, bottom=0.19, top=0.84,
                          hspace=0.54, wspace=0.32)
    axa, axb, axc, axd = [fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    bottom = np.zeros(len(prop))
    for i, col in enumerate(prop.columns):
        axa.bar(np.arange(len(prop)), prop[col], bottom=bottom, label=col.replace("_", " "), color=COLORS[i % len(COLORS)], width=0.78)
        bottom += prop[col].to_numpy()
    axa.axvline(4.5, color="#C7CDD3", lw=1.2, ls="--")
    axa.set_xticks(np.arange(len(prop)), [short_sample(x) for x in prop.index], rotation=45, ha="right")
    axa.set_ylabel("Cell composition (%)")
    axa.set_ylim(0, 100)
    axa.set_title("Conservative annotation by donor", loc="left")
    axa.legend(ncol=5, bbox_to_anchor=(0.0, -0.25), loc="upper left",
               fontsize=6.0, columnspacing=0.9, handlelength=1.6)
    label(axa, "A")

    group = counts.groupby(["condition", "conservative_label"], as_index=False)["n_cells"].sum()
    group["fraction"] = group["n_cells"] / group.groupby("condition")["n_cells"].transform("sum") * 100
    selected = list(label_totals.head(12).index)
    mat = group[group["conservative_label"].isin(selected)].pivot(index="conservative_label", columns="condition", values="fraction").fillna(0).reindex(selected)
    y = np.arange(len(mat))
    axb.barh(y + 0.18, mat.get("PBC", pd.Series(0, index=mat.index)), height=0.34, color=VERMILLION, label="PBC")
    axb.barh(y - 0.18, mat.get("control", pd.Series(0, index=mat.index)), height=0.34, color=BLUE, label="Control")
    axb.set_yticks(y, [x.replace("_", " ") for x in mat.index])
    axb.invert_yaxis(); axb.set_xlabel("Fraction of group cells (%)")
    axb.set_title("Group-level cell composition", loc="left"); axb.legend(loc="lower right")
    axb.text(0.99, 1.01, "Descriptive only; no differential-abundance inference",
             transform=axb.transAxes, ha="right", va="bottom", fontsize=6.5,
             color="#687681")
    style(axb, xgrid=True); label(axb, "B")

    chosen = [x for x in ["lrNK", "cNK", "CD4_T", "CD8_T", "TNK_mixed", "B_cell", "Myeloid", "Cholangiocyte"] if x in meta["conservative_label"].unique()]
    vals = [meta.loc[meta["conservative_label"] == c, "annotation_delta"].dropna().to_numpy() for c in chosen]
    bp = axc.boxplot(vals, orientation="horizontal", tick_labels=[x.replace("_", " ") for x in chosen], patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], COLORS):
        patch.set_facecolor(color); patch.set_alpha(0.65)
    axc.set_xlabel("SingleR annotation delta")
    axc.set_title("Annotation confidence by conservative label", loc="left")
    style(axc, xgrid=True); label(axc, "C")

    if cross.empty:
        axd.text(0.5, 0.5, "No T/NK refinement records", ha="center", va="center")
    else:
        im = axd.imshow(cross.to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=100)
        axd.set_xticks(np.arange(cross.shape[1]), [str(x).replace("_", " ") for x in cross.columns], rotation=45, ha="right")
        axd.set_yticks(np.arange(cross.shape[0]), [str(x).replace("_", " ") for x in cross.index])
        for i in range(cross.shape[0]):
            for j in range(cross.shape[1]):
                v = cross.iloc[i, j]
                if v >= 1:
                    axd.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.1, color="white" if v > 55 else INK)
        fig.colorbar(im, ax=axd, fraction=0.035, pad=0.02, label="Row fraction (%)")
    axd.set_xlabel("Pure T/NK sensitivity label", labelpad=7); axd.set_ylabel("Primary conservative label")
    axd.set_title("Primary versus sensitivity T/NK labels", loc="left")
    label(axd, "D")
    fig.suptitle("Supplementary Figure S7 | Target-independent liver cell annotation",
                 fontsize=14, fontweight="bold", x=0.08, y=0.975, ha="left", color=NAVY)
    fig.text(0.08, 0.938,
             "Conservative reference mapping and a target-excluded T/NK sensitivity classification.",
             fontsize=8.5, color="#687681")
    fig.text(0.54, 0.018,
             "Annotation guardrail: IL12RB2 and every gene in the prespecified candidate chain were excluded from annotation features.",
             fontsize=6.9, color=GREEN, ha="center", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F3FAF7", edgecolor="#CDE8DC"))
    save(fig, "FigS07_HRA008003_annotation_final")


def load_specs():
    rows = []
    pattern = re.compile(r"(?P<scheme>primary_conservative|pure_TNK_sensitivity)__(?P<filter>all_final_singlets|at_least_500_UMI)__(?P<cell>.+)_min(?P<minimum>10|20)_edgeR.tsv")
    for p in sorted((PB / "all_gene_results").glob("*.tsv")):
        m = pattern.fullmatch(p.name)
        if not m:
            continue
        d = read_tsv(p, usecols=lambda c: c in {"gene_symbol", "logFC", "PValue", "FDR"})
        hit = d[d["gene_symbol"] == "IL12RB2"]
        if hit.empty:
            continue
        row = m.groupdict() | hit.iloc[0].to_dict()
        rows.append(row)
    d = pd.DataFrame(rows)
    d["minimum"] = d["minimum"].astype(int)
    return d


def fig_s3():
    specs = load_specs()
    donor = read_tsv(PB / "04_donor_target_logCPM.tsv")
    donor = donor[(donor["gene_symbol"] == "IL12RB2") & (donor["scheme"] == "primary_conservative") & (donor["cell_filter"] == "all_final_singlets") & (donor["minimum_cells"] == 20)]
    perm = read_tsv(PB / "05_exact_donor_permutation_results.tsv")
    perm = perm[(perm["gene_symbol"] == "IL12RB2") & (perm["scheme"] == "primary_conservative") & (perm["cell_filter"] == "all_final_singlets") & (perm["minimum_cells"] == 20)]
    loo = read_tsv(PB / "07_target_leave_one_donor_out.tsv")
    loo = loo[loo["gene_symbol"] == "IL12RB2"]
    specs.to_csv(SRC / "FigS8_IL12RB2_specifications.tsv", sep="\t", index=False)
    donor.to_csv(SRC / "FigS8_IL12RB2_donor_logCPM.tsv", sep="\t", index=False)
    perm.to_csv(SRC / "FigS8_IL12RB2_exact_tests.tsv", sep="\t", index=False)
    loo.to_csv(SRC / "FigS8_IL12RB2_LODO.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.5, 9.8), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.98, bottom=0.12, top=0.84,
                          hspace=0.56, wspace=0.34)
    axa, axb, axc, axd = [fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    cell_order = [c for c in ["CD4_T", "cNK", "lrNK", "CD8_T", "TNK_mixed", "T_other"] if c in specs["cell"].unique()]
    spec_order = [("primary_conservative", "all_final_singlets", 20), ("primary_conservative", "all_final_singlets", 10),
                  ("primary_conservative", "at_least_500_UMI", 20), ("pure_TNK_sensitivity", "all_final_singlets", 20),
                  ("pure_TNK_sensitivity", "at_least_500_UMI", 20)]
    spec_names = ["Primary\nall · min20", "Primary\nall · min10", "Primary\n≥500 UMI · min20",
                  "Pure T/NK\nall · min20", "Pure T/NK\n≥500 UMI · min20"]
    mat = np.full((len(cell_order), len(spec_order)), np.nan)
    for i, cell in enumerate(cell_order):
        for j, (scheme, flt, minimum) in enumerate(spec_order):
            x = specs[(specs["cell"] == cell) & (specs["scheme"] == scheme) & (specs["filter"] == flt) & (specs["minimum"] == minimum)]
            if len(x): mat[i, j] = x.iloc[0]["logFC"]
    lim = max(2, np.nanmax(np.abs(mat)))
    im = axa.imshow(mat, aspect="auto", cmap=LinearSegmentedColormap.from_list("effect", [BLUE, "white", VERMILLION]), vmin=-lim, vmax=lim)
    axa.set_xticks(np.arange(len(spec_order)), spec_names, rotation=0, ha="center", fontsize=6.3)
    axa.set_yticks(np.arange(len(cell_order)), [x.replace("_", " ") for x in cell_order])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]): axa.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", fontsize=6.5)
    fig.colorbar(im, ax=axa, fraction=0.035, pad=0.02, label="PBC vs control log2FC")
    axa.set_title("Prespecified analysis-specification stability", loc="left")
    axa.text(0.99, 1.01, f"40/40 effects positive; minimum FDR={specs['FDR'].min():.3f}",
             transform=axa.transAxes, ha="right", va="bottom", fontsize=6.5,
             color=VERMILLION, fontweight="bold")
    label(axa, "A")

    positions, tickpos, ticklab = [], [], []
    cell_order2 = [c for c in ["CD4_T", "cNK", "lrNK", "CD8_T", "TNK_mixed"] if c in donor["cell_type"].unique()]
    for i, cell in enumerate(cell_order2):
        sub = donor[donor["cell_type"] == cell]
        for cond, off, color in [("control", -0.13, BLUE), ("PBC", 0.13, VERMILLION)]:
            vals = sub.loc[sub["condition"] == cond, "logCPM"].to_numpy()
            jitter = np.linspace(-0.045, 0.045, len(vals)) if len(vals) > 1 else np.array([0])
            axb.scatter(np.full(len(vals), i + off) + jitter, vals, color=color, s=28, edgecolor="white", lw=0.5, zorder=3)
            if len(vals): axb.plot([i + off - 0.09, i + off + 0.09], [np.mean(vals)]*2, color=color, lw=2)
        positions.append(i); tickpos.append(i); ticklab.append(cell.replace("_", " "))
    axb.set_xticks(tickpos, ticklab, rotation=25, ha="right")
    axb.set_ylabel("Donor IL12RB2 logCPM")
    axb.set_title("Donor-level expression, not cell-level pseudoreplication", loc="left")
    axb.legend(handles=[Line2D([0],[0], marker='o', color='none', markerfacecolor=BLUE, label='Control'), Line2D([0],[0], marker='o', color='none', markerfacecolor=VERMILLION, label='PBC')], loc="upper right")
    axb.text(0.01, 0.02, "Independent units: 5 PBC and 5 control donors",
             transform=axb.transAxes, fontsize=6.5, color="#687681")
    style(axb, ygrid=True); label(axb, "B")

    metrics = [x for x in ["logCPM", "log2_CPM_plus_0.5", "positive_fraction"] if x in perm["metric"].unique()]
    pm = perm.pivot(index="cell_type", columns="metric", values="exact_PValue").reindex(cell_order2)[metrics]
    im2 = axc.imshow(pm.to_numpy(), aspect="auto", cmap="YlOrRd_r", vmin=0, vmax=1)
    metric_names={"logCPM":"logCPM", "log2_CPM_plus_0.5":"log2(CPM + 0.5)",
                  "positive_fraction":"Positive-cell fraction"}
    axc.set_xticks(np.arange(len(metrics)), [metric_names.get(x,x) for x in metrics], rotation=0, ha="center", fontsize=6.7)
    axc.set_yticks(np.arange(len(pm)), [x.replace("_", " ") for x in pm.index])
    for i in range(pm.shape[0]):
        for j in range(pm.shape[1]):
            if np.isfinite(pm.iloc[i, j]): axc.text(j, i, f"{pm.iloc[i,j]:.3f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im2, ax=axc, fraction=0.035, pad=0.02, label="Exact two-sided P")
    axc.set_title("Exact enumeration across 252 case-label assignments", loc="left", pad=8)
    axc.text(0.99, 1.01, f"Minimum exact two-sided P={np.nanmin(pm.to_numpy()):.4f}; none <0.05",
             transform=axc.transAxes, ha="right", va="bottom", fontsize=6.5,
             color=VERMILLION, fontweight="bold")
    label(axc, "C")

    summary = loo.groupby("cell_type", as_index=False).agg(full=("full_effect_PBC_minus_control_logCPM", "first"), low=("leave_one_out_effect_logCPM", "min"), high=("leave_one_out_effect_logCPM", "max"), preserved=("same_direction_as_full", "mean"))
    summary = summary[summary["cell_type"].isin(cell_order2)].set_index("cell_type").reindex(cell_order2).reset_index()
    y = np.arange(len(summary))[::-1]
    axd.errorbar(summary["full"], y, xerr=[summary["full"]-summary["low"], summary["high"]-summary["full"]], fmt="o", color=GREEN, capsize=3)
    axd.axvline(0, color="#666", lw=0.9)
    axd.set_yticks(y, [x.replace("_", " ") for x in summary["cell_type"]])
    axd.set_xlabel("PBC − control donor mean logCPM")
    axd.set_title("Leave-one-donor-out effect ranges", loc="left")
    for yi, (_, r) in zip(y, summary.iterrows()): axd.text(r["high"] + 0.05, yi, f"{r['preserved']*100:.0f}% same direction", va="center", fontsize=6.5)
    style(axd, xgrid=True); label(axd, "D")
    axd.text(0.99, 0.02, "Robustness diagnostic—not independent replication",
             transform=axd.transAxes, ha="right", va="bottom", fontsize=6.5,
             color="#687681")
    fig.suptitle("Supplementary Figure S8 | HRA008003 donor-level pseudobulk sensitivity",
                 fontsize=14, fontweight="bold", x=0.08, y=0.975, ha="left", color=NAVY)
    fig.text(0.08, 0.938,
             "Directionally stable IL12RB2 effects remain non-confirmatory after multiplicity and exact donor permutation.",
             fontsize=8.5, color="#687681")
    fig.text(0.54, 0.026,
             "Inference boundary: no IL12RB2 contrast passed transcriptome-wide FDR 0.05; sensitivity analyses do not create independent replication.",
             fontsize=7.0, color=VERMILLION, ha="center", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF7F2", edgecolor="#F3D5C4"))
    save(fig, "FigS08_HRA008003_pseudobulk_sensitivity_final")


def fig_s4():
    bulk = read_tsv(BULK / "16_IL12RB2_all_cohort_comparisons.tsv")
    effects = read_tsv(PATHWAY / "03_all_cohort_effects.tsv")
    meta = read_tsv(PATHWAY / "04_primary_random_effect_meta.tsv")
    lodo = read_tsv(PATHWAY / "05_leave_one_dataset_out.tsv")
    logo = read_tsv(PATHWAY / "06_leave_one_gene_out.tsv")
    trim = read_tsv(PATHWAY / "07_inflammatory_driver_trimmed_meta.tsv")
    bulk.to_csv(SRC / "FigS9_bulk_IL12RB2.tsv", sep="\t", index=False)
    effects.to_csv(SRC / "FigS9_pathway_cohort_effects.tsv", sep="\t", index=False)
    pd.concat([meta.assign(sensitivity="primary"), lodo.assign(sensitivity="LODO"), logo.assign(sensitivity="LOGO"), trim.assign(sensitivity="trimmed")], ignore_index=True, sort=False).to_csv(SRC / "FigS9_pathway_sensitivity.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.5, 9.6), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.98, bottom=0.12, top=0.84,
                          hspace=0.45, wspace=0.42)
    axa, axb, axc, axd = [fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    b = bulk.copy(); y = np.arange(len(b))[::-1]
    for yi, (_, r) in zip(y, b.iterrows()):
        color = VERMILLION if r["contrast"] in ["PBC_vs_control", "treatment_naive_vs_control"] else GRAY
        low, high = r["hedges_g"] - 1.96*r["hedges_g_se"], r["hedges_g"] + 1.96*r["hedges_g_se"]
        axa.errorbar(r["hedges_g"], yi, xerr=[[r["hedges_g"]-low],[high-r["hedges_g"]]], fmt="o", color=color, capsize=3)
    axa.axvline(0, color="#666", lw=0.9)
    labels = [f"{r.dataset.replace('_',' ')}\n{r.contrast.replace('_',' ')}" for r in b.itertuples()]
    axa.set_yticks(y, labels); axa.set_xlabel("IL12RB2 Hedges g (95% CI)")
    axa.set_title("Bulk IL12RB2 effects across available contrasts", loc="left")
    axa.legend(handles=[Line2D([0],[0], marker="o", color="none", markerfacecolor=VERMILLION,
                               label="Primary PBC/disease-context contrasts"),
                         Line2D([0],[0], marker="o", color="none", markerfacecolor=GRAY,
                               label="Secondary clinical/disease context")],
               loc="upper left", bbox_to_anchor=(0.0, -0.18), ncol=2,
               fontsize=6.2, handletextpad=0.3, columnspacing=1.0)
    primary_bulk = b[b["contrast"].isin(["PBC_vs_control", "treatment_naive_vs_control"])]
    available_fdr = primary_bulk["limma_fdr"].notna().sum()
    passed_fdr = (primary_bulk["limma_fdr"].dropna() < 0.05).sum()
    axa.text(0.99, 1.01,
             f"Primary contrasts: {int((primary_bulk['hedges_g']>0).sum())}/{len(primary_bulk)} positive; {passed_fdr}/{available_fdr} available FDR<0.05",
             transform=axa.transAxes, ha="right", va="bottom", fontsize=6.3,
             color=VERMILLION, fontweight="bold")
    style(axa, xgrid=True); label(axa, "A")

    def pathway_forest(ax, signature, letter):
        e = effects[(effects["signature"] == signature) & (effects["score_method"] == "robust_z") & effects["contrast"].str.contains("PBC_vs_healthy|treatment_naive_vs_control", regex=True)].copy()
        m = meta[(meta["signature"] == signature) & (meta["score_method"] == "robust_z")].iloc[0]
        yy = np.arange(len(e)+1)[::-1]
        for yi, (_, r) in zip(yy[:-1], e.iterrows()):
            ax.errorbar(r["hedges_g"], yi, xerr=1.96*r["se"], fmt="o", color=BLUE, capsize=3)
        ax.errorbar(m["estimate"], yy[-1], xerr=[[m["estimate"]-m["ci_low_mKH"]],[m["ci_high_mKH"]-m["estimate"]]], fmt="D", color=VERMILLION, capsize=3, ms=5.5)
        ax.axvline(0, color="#666", lw=0.9)
        ax.set_yticks(yy, [x.replace("_", " ") for x in e["dataset"]] + ["mKH random-effects meta"])
        ax.set_xlabel("Hedges g (95% CI)")
        ax.set_title(signature.replace("_", "–") + " program", loc="left")
        ax.text(0.98, 0.05,
                f"4/4 positive | pooled P={m['p_mKH']:.3f}; I²={m['I2_percent']:.1f}%",
                transform=ax.transAxes, ha="right", fontsize=6.7, color=VERMILLION,
                fontweight="bold")
        style(ax, xgrid=True); label(ax, letter)
    pathway_forest(axb, "IL12_STAT4", "B")
    pathway_forest(axc, "IFNG_RESPONSE", "C")

    rows = []
    for sig in ["IL12_STAT4", "IFNG_RESPONSE"]:
        mm = meta[(meta["signature"] == sig) & (meta["score_method"] == "robust_z")].iloc[0]
        rows.append((sig, "Primary", mm["estimate"], mm["ci_low_mKH"], mm["ci_high_mKH"]))
        for _, r in lodo[lodo["signature"] == sig].iterrows(): rows.append((sig, "omit " + r["left_out_dataset"].replace("_", " "), r["estimate"], r["ci_low_mKH"], r["ci_high_mKH"]))
        tr = trim[trim["signature"].str.startswith(sig)]
        if len(tr):
            r = tr.iloc[0]; rows.append((sig, "trimmed genes", r["estimate"], r["ci_low_mKH"], r["ci_high_mKH"]))
    sd = pd.DataFrame(rows, columns=["signature","analysis","estimate","low","high"])
    y = np.arange(len(sd))[::-1]
    for yi, (_, r) in zip(y, sd.iterrows()):
        color = VERMILLION if r["analysis"] == "Primary" else (BLUE if r["signature"] == "IL12_STAT4" else GREEN)
        axd.errorbar(r["estimate"], yi, xerr=[[r["estimate"]-r["low"]],[r["high"]-r["estimate"]]], fmt="o", color=color, capsize=2.5, ms=4)
    axd.axvline(0, color="#666", lw=0.9)
    sig_short={"IL12_STAT4":"IL12–STAT4","IFNG_RESPONSE":"IFNG"}
    clean_analysis=[]
    for r in sd.itertuples():
        a=str(r.analysis).replace("GSE119600 whole blood","GSE119600").replace("GSE159676 liver","GSE159676")
        a=a.replace("GSE79850 liver NanoString","GSE79850").replace("GSE138888 DNT","GSE138888")
        clean_analysis.append(f"{sig_short.get(r.signature,r.signature)}: {a}")
    axd.set_yticks(y, clean_analysis, fontsize=6.8)
    axd.set_xlabel("Sensitivity meta estimate (95% mKH CI)")
    axd.set_title("Dataset and signature sensitivity analyses", loc="left")
    axd.text(0.99, 0.02, "All displayed 95% mKH intervals cross zero",
             transform=axd.transAxes, ha="right", va="bottom", fontsize=6.5,
             color="#687681")
    style(axd, xgrid=True); label(axd, "D")
    fig.suptitle("Supplementary Figure S9 | Bulk transcriptomic and pathway context",
                 fontsize=14, fontweight="bold", x=0.10, y=0.975, ha="left", color=NAVY)
    fig.text(0.10, 0.938,
             "Cross-cohort directional consistency is separated from confirmatory significance and heterogeneity.",
             fontsize=8.5, color="#687681")
    fig.text(0.55, 0.025,
             "Inference boundary: bulk data support directional context, not an IL12RB2 diagnostic/prognostic biomarker; pathway meta-analyses remain non-confirmatory.",
             fontsize=7.0, color=VERMILLION, ha="center", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF7F2", edgecolor="#F3D5C4"))
    save(fig, "FigS09_bulk_pathway_meta_final")


def fig_s5():
    robust = read_tsv(PXD / "11_PXD083842_TYK2_robust_methods.tsv")
    robust = robust[robust["contrast"] == "PBC_vs_Healthy"].copy()
    loo = read_tsv(PXD / "10_PXD083842_TYK2_leave_one_out.tsv")
    loo = loo[loo["contrast"] == "PBC_vs_Healthy"].copy()
    audit = read_tsv(PXD / "17_PXD083842_covariate_metadata_audit.tsv")
    adjusted = read_tsv(PXD / "18_PXD083842_TYK2_adjusted_results.tsv")
    robust.to_csv(SRC / "FigS10_TYK2_robust_methods.tsv", sep="\t", index=False)
    loo.to_csv(SRC / "FigS10_TYK2_LOSO.tsv", sep="\t", index=False)
    audit.to_csv(SRC / "FigS10_confounder_audit.tsv", sep="\t", index=False)
    adjusted.to_csv(SRC / "FigS10_adjusted_evidence.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.5, 9.7), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.98, bottom=0.13, top=0.84,
                          hspace=0.45, wspace=0.38)
    axa, axb, axc, axd = [fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    a = adjusted.copy(); y = np.arange(len(a))[::-1]
    for yi, (_, r) in zip(y, a.iterrows()):
        if pd.notna(r["log2FC"]):
            color = VERMILLION if bool(r["TYK2_in_FDR_lt_0_05_list"]) else GRAY
            axa.scatter(r["log2FC"], yi, color=color, s=52, zorder=3)
            axa.text(r["log2FC"]+0.015, yi, f"FDR={r['adjusted_p']:.3f}", va="center", fontsize=6.5)
        else:
            axa.scatter(0, yi, marker="x", color=GRAY, s=44)
            axa.text(0.015, yi, "not in released FDR<0.05 list; β unavailable", va="center", fontsize=6.2, color="#666")
    axa.axvline(0, color="#666", lw=0.9)
    axa.set_yticks(y, [x.replace("_", " ") for x in a["adjustment_or_exclusion"]])
    axa.set_xlabel("Reported TYK2 log2FC where available")
    axa.set_title("Adjustment and exclusion audit", loc="left")
    axa.text(0.99, 1.01, "Only 2/8 exact TYK2 coefficients were publicly released",
             transform=axa.transAxes, ha="right", va="bottom", fontsize=6.5,
             color=VERMILLION, fontweight="bold")
    style(axa, xgrid=True); label(axa, "A")

    method_order = ["Welch","Mann_Whitney","Permutation_100000","Hodges_Lehmann","OLS_HC3","Huber_RLM","winsor_5_95","MAD_3sigma"]
    r = robust.set_index("method").reindex(method_order).dropna(how="all").reset_index(); y = np.arange(len(r))[::-1]
    axb.scatter(r["estimate"], y, c=[VERMILLION if p < 0.05 else GRAY for p in r["p_value"]], s=48)
    for yi, (_, z) in zip(y, r.iterrows()): axb.text(z["estimate"]+0.006, yi, f"P={z['p_value']:.3g}", va="center", fontsize=6.5)
    axb.axvline(0, color="#666", lw=0.9)
    axb.set_yticks(y, [x.replace("_", " ") for x in r["method"]])
    axb.set_xlabel("PBC − healthy TYK2 abundance")
    axb.set_title("Distributional and model robustness", loc="left")
    axb.text(0.99, 1.01,
             f"8/8 positive; estimate {r['estimate'].min():.3f}–{r['estimate'].max():.3f}",
             transform=axb.transAxes, ha="right", va="bottom", fontsize=6.5,
             color=VERMILLION, fontweight="bold")
    style(axb, xgrid=True); label(axb, "B")

    colors = [VERMILLION if x == "PBC" else BLUE for x in loo["dropped_group"]]
    axc.scatter(np.arange(len(loo)), loo["delta"], c=colors, s=33, edgecolor="white", lw=0.4)
    axc.axhline(robust.loc[robust["method"] == "Welch", "estimate"].iloc[0], color=GREEN, ls="--", lw=1.2, label="Full-sample estimate")
    axc.axhline(0, color="#666", lw=0.8)
    axc.set_xticks(np.arange(len(loo)), loo["dropped_sample"], rotation=70, fontsize=6)
    axc.set_ylabel("TYK2 PBC − healthy difference")
    axc.set_title("Leave-one-sample-out stability", loc="left")
    axc.legend(handles=[Line2D([0],[0], color=GREEN, ls="--", lw=1.2, label="Full-sample estimate"),
                        Line2D([0],[0], marker="o", color="none", markerfacecolor=VERMILLION,
                               label="Dropped PBC sample"),
                        Line2D([0],[0], marker="o", color="none", markerfacecolor=BLUE,
                               label="Dropped healthy sample")],
               loc="upper left", fontsize=6.2)
    axc.text(0.99, 1.01, f"{len(loo)}/{len(loo)} deletions retained a positive effect",
             transform=axc.transAxes, ha="right", va="bottom", fontsize=6.5,
             color=VERMILLION, fontweight="bold")
    style(axc, ygrid=True); label(axc, "C")

    display = audit.copy()
    display["risk_score"] = display["risk"].map({"critical": 3, "high": 2, "moderate": 1, "low": 0, "low_for_TYK2": 1}).fillna(1)
    display["patient_public"] = display["patient_level_public"].astype(str).str.lower().eq("true").astype(int)
    support_text = display["TYK2_result_after_control"].fillna("").str.lower()
    display["TYK2_FDR_support"] = (support_text.str.contains(r"fdr\s*0\.0162|remained in|retained in", regex=True) & ~support_text.str.contains(r"not fdr|not present|no liver", regex=True)).astype(int)
    matrix = display[["risk_score", "patient_public", "TYK2_FDR_support"]].to_numpy()
    rgb = np.ones((len(display), 3, 3))
    risk_palette={0:"#E8F5E9",1:"#F0E442",2:"#E69F00",3:"#D55E00"}
    from matplotlib.colors import to_rgb
    for i in range(len(display)):
        rgb[i,0,:]=to_rgb(risk_palette[int(matrix[i,0])])
        rgb[i,1,:]=to_rgb("#D9EEF8" if matrix[i,1] else "#F0F2F4")
        rgb[i,2,:]=to_rgb("#D9F0E5" if matrix[i,2] else "#F0F2F4")
    axd.imshow(rgb, aspect="auto")
    axd.set_yticks(np.arange(len(display)), [str(x).replace("_", " ") for x in display["variable"]])
    axd.set_xticks([0,1,2], ["Confounding risk\n0–3", "Covariate publicly\nobservable", "TYK2 retained in\nreleased FDR list"])
    for i in range(matrix.shape[0]):
        axd.text(0, i, f"{int(matrix[i,0])}", ha="center", va="center", fontsize=7)
        axd.text(1, i, "Yes" if matrix[i,1] else "No", ha="center", va="center", fontsize=7)
        axd.text(2, i, "Yes" if matrix[i,2] else "No", ha="center", va="center", fontsize=7)
    axd.set_title("Structural confounding and identification limits", loc="left")
    label(axd, "D")
    fig.suptitle("Supplementary Figure S10 | PXD083842 proteomic robustness and confounding",
                 fontsize=14, fontweight="bold", x=0.10, y=0.975, ha="left", color=NAVY)
    fig.text(0.10, 0.938,
             "Positive TYK2 abundance is distributionally robust, but disease specificity is limited by structural confounding.",
             fontsize=8.5, color="#687681")
    fig.text(0.55, 0.026,
             "Inference boundary: age and recruitment source are structurally confounded with disease; unavailable coefficients and FDR-list absence cannot be interpreted as effect reversal.",
             fontsize=7.0, color=VERMILLION, ha="center", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF7F2", edgecolor="#F3D5C4"))
    save(fig, "FigS10_PXD083842_robustness_confounding_final")


def fig_s6():
    summary = read_tsv(TYK2 / "07_TYK2_LD_bootstrap_summary.tsv")
    pop = read_tsv(TYK2 / "06_TYK2_population_LD_sensitivity.tsv")
    diff = read_tsv(TYK2 / "08b_TYK2_FIN_vs_nonFIN_LD_difference.tsv")
    draws_path = TYK2 / "08_TYK2_LD_bootstrap_draws.tsv.gz"
    draws = read_tsv(draws_path, compression="gzip")
    valid = draws[draws["valid"] == True].copy()
    summary.to_csv(SRC / "FigS11_LD_bootstrap_summary.tsv", sep="\t", index=False)
    pop.to_csv(SRC / "FigS11_population_sensitivity.tsv", sep="\t", index=False)
    diff.to_csv(SRC / "FigS11_FIN_nonFIN_difference.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.5, 9.6))
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.97, bottom=0.13, top=0.86,
                          hspace=0.42, wspace=0.34)
    axa, axb, axc, axd = [fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    scenarios = list(summary["scenario"])
    names = {"EUR_all_unstratified":"EUR unstratified", "EUR_all_population_stratified":"EUR stratified", "FIN":"FIN", "EUR_non_FIN_population_stratified":"EUR non-FIN"}
    for i, sc in enumerate(scenarios):
        v = valid.loc[valid["scenario"] == sc, "signed_r"].to_numpy()
        axa.hist(v, bins=70, density=True, histtype="step", lw=1.5, color=COLORS[i], label=names.get(sc, sc))
    axa.set_xlabel("Bootstrap signed dosage correlation r")
    axa.set_ylabel("Density")
    axa.set_title("Reference-LD bootstrap distributions", loc="left")
    axa.legend(ncol=2, loc="upper right", handlelength=2.2, columnspacing=1.0)
    style(axa, xgrid=True); label(axa, "A")

    centers = np.arange(len(scenarios)) + 1
    offsets = [-0.18, 0.18]
    variant_specs = [
        ("rs345_conditional_p", "rs34536443", VERMILLION),
        ("rs230_conditional_p", "rs2304256", BLUE),
    ]
    for offset, (column, variant, color) in zip(offsets, variant_specs):
        data = [
            -np.log10(valid.loc[valid["scenario"] == sc, column].clip(lower=1e-300))
            for sc in scenarios
        ]
        parts = axb.violinplot(data, positions=centers + offset, widths=0.31,
                               showmeans=False, showmedians=True, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(0.68)
        parts["cmedians"].set_color(INK)
        parts["cmedians"].set_linewidth(1.1)
    axb.axhline(-np.log10(0.05), color="#666", ls="--", lw=0.9)
    axb.text(4.40, -np.log10(0.05)+0.10, "P=0.05", ha="right", va="bottom",
             fontsize=6.7, color="#555")
    axb.set_xticks(centers, ["EUR\nunstratified", "EUR\nstratified", "FIN", "EUR\nnon-FIN"],
                   rotation=0, ha="center", fontsize=7.2)
    axb.set_xlim(0.55, 4.45)
    axb.set_ylabel("Conditional −log10(P)")
    axb.set_title("Both signals persist across LD uncertainty", loc="left")
    axb.legend(handles=[
        Line2D([0], [0], marker="s", color="none", markerfacecolor=VERMILLION,
               markeredgecolor=VERMILLION, markersize=7, label="rs34536443"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=BLUE,
               markeredgecolor=BLUE, markersize=7, label="rs2304256"),
    ], ncol=2, loc="upper right", handletextpad=0.4, columnspacing=0.9)
    both = [float(summary.loc[summary["scenario"] == sc, "fraction_both_conditional_p_lt_0_05"].iloc[0]) for sc in scenarios]
    axb.text(0.02, 0.04,
             "Both conditional P<0.05: " + " · ".join(f"{x*100:.3f}%" if x < 0.99999 else "100%" for x in both),
             transform=axb.transAxes, fontsize=6.6, color=INK,
             bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor=LIGHTGRAY))
    style(axb, ygrid=True); label(axb, "B")

    s = summary.copy(); y = np.arange(len(s))[::-1]
    axc.errorbar(s["observed_r"], y, xerr=[s["observed_r"]-s["r_CI95_low"], s["r_CI95_high"]-s["observed_r"]], fmt="o", color=GREEN, capsize=3)
    axc.axvline(0, color="#666", lw=0.9)
    axc.set_yticks(y, [names.get(x,x) for x in s["scenario"]])
    axc.set_xlabel("Observed signed r and bootstrap 95% CI")
    axc.set_title("Population and resampling sensitivity", loc="left")
    for yi, (_, r) in zip(y, s.iterrows()):
        suffix = " | 3 C alleles" if r["scenario"] == "FIN" else ""
        axc.text(r["r_CI95_high"]+0.012, yi,
                 f"r={r['observed_r']:.3f} | valid={r['valid_fraction']*100:.1f}%{suffix}",
                 va="center", fontsize=6.3)
    axc.set_xlim(-0.02, 0.52)
    style(axc, xgrid=True); label(axc, "C")

    invalid = (1-summary["valid_fraction"])*100
    bars = axd.bar(np.arange(len(summary)), invalid,
                   color=["#C8D0D8", "#C8D0D8", GREEN, "#C8D0D8"], width=0.58)
    axd.set_xticks(np.arange(len(summary)), ["EUR\nunstratified", "EUR\nstratified", "FIN", "EUR\nnon-FIN"])
    axd.set_ylabel("Invalid monomorphic replicates (%)")
    axd.set_title("Rare-allele instability is concentrated in FIN", loc="left")
    for bar, value in zip(bars, invalid):
        axd.text(bar.get_x()+bar.get_width()/2, max(value, 0)+0.07,
                 f"{value:.3f}%" if value else "0%", ha="center", va="bottom", fontsize=6.8)
    d = diff.iloc[0]
    axd.text(0.02, 0.95,
             f"FIN − non-FIN Δr = {d['observed_r_difference']:.3f}\n"
             f"95% CI {d['difference_CI95_low']:.3f} to {d['difference_CI95_high']:.3f}\n"
             f"bootstrap P = {d['two_sided_bootstrap_p_around_zero']:.3f}\n"
             "No significant population difference",
             transform=axd.transAxes, ha="left", va="top", fontsize=7.3,
             bbox=dict(facecolor="white", edgecolor=LIGHTGRAY, boxstyle="round,pad=0.35"))
    axd.text(0.02, 0.58, "FIN: n=99; 3 rs34536443-C alleles\n4,767/100,000 replicates monomorphic",
             transform=axd.transAxes, ha="left", va="top", fontsize=6.8, color=INK,
             bbox=dict(facecolor="#F3FBF8", edgecolor="#B9E1D3", boxstyle="round,pad=0.30"))
    style(axd, ygrid=True); label(axd, "D")
    fig.suptitle("Supplementary Figure S11 | TYK2 LD bootstrap and FIN/non-FIN sensitivity",
                 fontsize=14, fontweight="bold", x=0.10, y=0.965, ha="left", color=NAVY)
    fig.text(0.10, 0.928,
             "Reference-LD resampling supports two-signal separation while exposing rare-allele instability in FIN.",
             fontsize=8.5, color="#687681")
    fig.text(0.535, 0.035,
             "Inference boundary: bootstrap propagates sampling uncertainty in 1000 Genomes reference LD; 1000G FIN is not FinnGen study-specific signed LD and does not replace individual-level conditional regression.",
             fontsize=6.9, color=VERMILLION, ha="center", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF7F2", edgecolor="#F3D5C4"))
    save(fig, "FigS11_TYK2_LD_bootstrap_final")


def main():
    fig_s1()
    fig_s2()
    fig_s3()
    fig_s4()
    fig_s5()
    fig_s6()
    print(f"Created renumbered Supplementary Figures S6-S11 in {OUT}")


if __name__ == "__main__":
    main()
