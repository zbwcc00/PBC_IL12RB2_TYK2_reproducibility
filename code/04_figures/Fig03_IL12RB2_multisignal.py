from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch
import numpy as np
from scipy.stats import norm


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "figures" / "main"
ANALYSIS = REPO / "data" / "derived" / "analysis"
DISCOVERY = ANALYSIS / "genetics" / "discovery"
REPLICATION = ANALYSIS / "genetics" / "replication"
DISC = DISCOVERY / "21_coloc_susie_signal_pairs_l3_review.tsv"
DISC_SENS = DISCOVERY / "22_coloc_susie_prior_sensitivity_l3_review.tsv"
EXT = REPLICATION / "03_external_signal_pair_primary.tsv"
FG_REP = REPLICATION / "05_finngen_prespecified_variant_replication.tsv"
FG_REGION = REPLICATION / "04_finngen_region_qc_top_signal.tsv"
FG_ABF = REPLICATION / "16_il12rb2_finngen_single_causal_abf_95cs.tsv"
BF_BF = REPLICATION / "17_il12rb2_finngen_external_eqtl_bf_bf_sensitivity.tsv"
DIRECTION = ANALYSIS / "genetics" / "direction" / "01_il12rb2_alt_aligned_effect_direction.tsv"
GENE = (
    ANALYSIS
    / "genetics"
    / "independent_replication"
    / "finngen_gene"
    / "FinnGen_R12_PBC_IL12RB2_gene_level.tsv"
)

BLUE, SKY, GREEN = "#0072B2", "#56B4E9", "#009E73"
ORANGE, RED, GRAY = "#E69F00", "#D55E00", "#A7B1B8"
INK, MUTED = "#263746", "#667784"


def read_tsv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def style_axis(ax, grid=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis=grid, color="#D8DFE4", lw=0.6, alpha=0.75)
    ax.set_axisbelow(True)


def panel_label(ax, label):
    ax.text(-0.10, 1.07, label, transform=ax.transAxes, ha="right", va="top",
            fontsize=14, fontweight="bold", color=INK)


def fmt_prob(value):
    if value >= 0.995:
        return f"{value:.3f}"
    if value >= 0.01:
        return f"{value:.2f}"
    return "<0.01"


def ci_from_beta_p(beta, pvalue):
    z = norm.isf(max(pvalue, 1e-300) / 2)
    se = abs(beta) / z
    return beta - 1.96 * se, beta + 1.96 * se


def rounded_box(ax, xy, width, height, face, edge, text):
    box = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.018,rounding_size=0.025",
                         facecolor=face, edgecolor=edge, linewidth=1.1,
                         transform=ax.transAxes, clip_on=False)
    ax.add_patch(box)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.25)


def main():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8.4,
        "axes.titlesize": 10.2, "axes.titleweight": "bold", "axes.labelsize": 9,
        "xtick.labelsize": 7.8, "ytick.labelsize": 7.8,
        "legend.fontsize": 7.1, "legend.frameon": False,
        "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(13.8, 9.2), facecolor="white")
    outer = fig.add_gridspec(2, 1, height_ratios=[0.96, 1.04], hspace=0.42)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.14], wspace=0.36)
    bottom = outer[1].subgridspec(1, 3, width_ratios=[1.00, 1.12, 1.08], wspace=0.64)
    agrid = top[0, 0].subgridspec(1, 2, width_ratios=[0.67, 0.33], wspace=0.38)
    ax_a = fig.add_subplot(agrid[0, 0])
    ax_a_prior = fig.add_subplot(agrid[0, 1])
    ax_b = fig.add_subplot(top[0, 1])
    ax_c, ax_d, ax_e = (fig.add_subplot(bottom[0, i]) for i in range(3))

    cmap = LinearSegmentedColormap.from_list("coloc", ["#F3F5F6", "#D9EAF1", SKY, BLUE, "#004B73"])

    # A | Discovery coloc.susie and prior sensitivity inset.
    disc = [r for r in read_tsv(DISC) if r["analysis_id"] == "SA09"]
    mat = np.zeros((3, 3))
    for row in disc:
        mat[int(row["pbc_signal"]) - 1, int(row["eqtl_signal"]) - 1] = float(row["PP.H4"])
    im = ax_a.imshow(mat, vmin=0, vmax=1, cmap=cmap, aspect="auto")
    for i in range(3):
        for j in range(3):
            val = mat[i, j]
            ax_a.text(j, i, fmt_prob(val), ha="center", va="center", fontsize=9,
                      fontweight="bold", color="white" if val > 0.55 else INK)
    ax_a.add_patch(plt.Rectangle((-0.48, -0.48), 0.96, 0.96, fill=False, ec=RED, lw=2))
    ax_a.set_xticks(range(3), ["eQTL S1", "eQTL S2", "eQTL S3"])
    ax_a.set_yticks(range(3), ["PBC S1", "PBC S2", "PBC S3"])
    ax_a.set_xlabel("IL12RB2/NK eQTL component")
    ax_a.set_ylabel("PBC component")
    ax_a.set_title("Discovery coloc.susie", loc="left", color=INK, pad=8)
    for spine in ax_a.spines.values(): spine.set_visible(False)
    ax_a.tick_params(length=0)
    sens = [r for r in read_tsv(DISC_SENS) if r["analysis_id"] == "SA09" and r["pbc_signal"] == "1" and r["eqtl_signal"] == "1"]
    sens.sort(key=lambda r: float(r["p12"]))
    sx = [float(r["p12"]) for r in sens]
    ax_a_prior.plot(sx, [float(r["PP.H4"]) for r in sens], color=BLUE, marker="o", ms=4.2, lw=1.6, label="H4")
    ax_a_prior.plot(sx, [float(r["PP.H3"]) for r in sens], color=RED, marker="s", ms=3.8, lw=1.4, label="H3")
    ax_a_prior.set_xscale("log"); ax_a_prior.set_ylim(0, 1.04)
    ax_a_prior.set_xticks([1e-6, 1e-5], ["1e-6", "1e-5"])
    ax_a_prior.set_xlabel("Shared-signal prior (p12)", fontsize=6.6)
    ax_a_prior.set_ylabel("Posterior", fontsize=6.6)
    ax_a_prior.set_title("Prior sensitivity", fontsize=7.2, loc="left", color=INK, pad=5)
    ax_a_prior.tick_params(labelsize=6.2)
    ax_a_prior.legend(fontsize=6.2, loc="center right")
    style_axis(ax_a_prior, "y")
    ax_a.text(0.01, -0.21, "Shared rs6679356: PP.H4=0.971–0.997; convergent L=3 eQTL review.",
              transform=ax_a.transAxes, fontsize=6.9, color=MUTED)
    panel_label(ax_a, "A")

    # B | External signal-level replication.
    ext = [r for r in read_tsv(EXT) if r["gene_symbol"] == "IL12RB2"]
    studies = [("Perez_2022", "Perez NK"), ("GTEx_v10", "GTEx blood"), ("Nathan_2022", "Nathan Th1")]
    emat, h3 = np.zeros((3, 3)), np.zeros((3, 3))
    for i, (study, _) in enumerate(studies):
        for row in [r for r in ext if r["study_label"] == study]:
            j = int(row["pbc_signal"]) - 1
            emat[i, j], h3[i, j] = float(row["PP.H4"]), float(row["PP.H3"])
    im2 = ax_b.imshow(emat, vmin=0, vmax=1, cmap=cmap, aspect="auto")
    for i in range(3):
        for j in range(3):
            val = emat[i, j]; fg = "white" if val > 0.55 else INK
            ax_b.text(j, i - 0.07, fmt_prob(val), ha="center", va="center", fontsize=9,
                      fontweight="bold", color=fg)
            ax_b.text(j, i + 0.22, f"H3 {h3[i,j]:.2f}", ha="center", va="center", fontsize=6.1,
                      color=fg if val > 0.55 else MUTED)
    for i in (0, 1):
        ax_b.add_patch(plt.Rectangle((-0.48, i - 0.48), 0.96, 0.96, fill=False, ec=GREEN, lw=2))
    ax_b.set_xticks(range(3), ["PBC S1", "PBC S2", "PBC S3"])
    ax_b.set_yticks(range(3), [x[1] for x in studies])
    ax_b.set_xlabel("PBC fine-mapped component")
    ax_b.set_ylabel("Independent molecular-QTL study")
    ax_b.set_title("External eQTL signal-level replication", loc="left", color=INK, pad=8)
    for spine in ax_b.spines.values(): spine.set_visible(False)
    ax_b.tick_params(length=0)
    cb2 = fig.colorbar(im2, ax=ax_b, fraction=0.038, pad=0.03); cb2.set_label("PP.H4", labelpad=6); cb2.set_ticks([0, .5, 1])
    ax_b.text(0.01, -0.21, "Perez NK and GTEx blood replicate S1; Nathan Th1 defines a context boundary.",
              transform=ax_b.transAxes, fontsize=6.9, color=MUTED)
    panel_label(ax_b, "B")

    # C | Disease effect forest plus FinnGen evidence ladder.
    direction = read_tsv(DIRECTION)
    disease = [r for r in direction if r["evidence_layer"].startswith("PBC_GWAS")]
    labels = ["Cordell", "FinnGen", "UK Biobank"]
    y = np.arange(3)[::-1]
    for yy, label, row, color in zip(y, labels, disease, [GRAY, BLUE, SKY]):
        beta, pval = float(row["beta_effect_allele"]), float(row["pvalue"])
        low, high = ci_from_beta_p(beta, pval)
        ax_c.errorbar(beta, yy, xerr=[[beta-low],[high-beta]], fmt="o", ms=5.5, capsize=2.5,
                      color=color, ecolor=color, lw=1.3)
        ax_c.text(-0.73, yy - .24, f"P={pval:.1e}", color=color, fontsize=6.3)
    ax_c.axvline(0, color="#6B737A", lw=.9); ax_c.set_xlim(-.77,.05); ax_c.set_ylim(-1.35,2.55)
    ax_c.set_yticks(y, labels); ax_c.set_xlabel("rs6679356-T effect on log-odds")
    ax_c.set_title("FinnGen replication at three resolutions", loc="left", color=INK, pad=8)
    style_axis(ax_c)
    gene = read_tsv(GENE)[0]
    rounded_box(ax_c,(0.00,-.33),.30,.18,"#E8F3F8",BLUE,"VARIANT\nP=2.55e-8")
    rounded_box(ax_c,(0.35,-.33),.30,.18,"#EAF6F1",GREEN,"REGION\nP=2.78e-9")
    rounded_box(ax_c,(0.70,-.33),.30,.18,"#FFF5E3",ORANGE,"GENE rank 91\nFDR=0.00484")
    ax_c.text(.70,-.40,"Bonferroni P=0.440",transform=ax_c.transAxes,fontsize=6.2,color=RED)
    panel_label(ax_c,"C")

    # D | Wakefield ABF + FinnGen-eQTL compatibility.
    abf = read_tsv(FG_ABF); priors=[.15,.20,.30]; sizes=[]; post=[]
    for prior in priors:
        rr=[r for r in abf if abs(float(r["prior_sd_key"])-prior)<1e-9]
        sizes.append(len(rr)); post.append(float(min(rr,key=lambda r:int(r["rank"]))["posterior"]))
    x=np.arange(3); ax_d2=ax_d.twinx()
    ax_d.bar(x,sizes,width=.55,color="#C9D5DC",edgecolor="white")
    ax_d2.plot(x,post,color=BLUE,marker="o",ms=5,lw=1.8)
    for xx,ss,pp in zip(x,sizes,post):
        ax_d.text(xx-.10,ss+.5,str(ss),ha="center",fontsize=6.6)
        ax_d2.text(xx+.10,pp+.010,f"{pp:.3f}",ha="center",fontsize=6.2,color=BLUE)
    ax_d.set_xticks(x,["0.15","0.20","0.30"]); ax_d.set_ylim(0,30); ax_d2.set_ylim(0,.16)
    ax_d.set_xlabel("Wakefield prior SD"); ax_d.set_ylabel("95% CS size")
    ax_d2.set_ylabel("Lead posterior", color=BLUE, rotation=270, labelpad=12, fontsize=6.8)
    ax_d2.set_yticks(np.arange(.02, .161, .02))
    ax_d2.tick_params(axis="y",colors=BLUE,labelsize=6.2,pad=1,length=0,right=False)
    ax_d.set_title("ABF and FinnGen–eQTL compatibility",loc="left",color=INK,pad=8)
    style_axis(ax_d); ax_d2.spines["top"].set_visible(False); ax_d2.spines["left"].set_visible(False); ax_d2.spines["right"].set_visible(False)
    ins=ax_d.inset_axes([.42,.17,.54,.48]); bfbf=read_tsv(BF_BF)
    inset_endpoints = []
    for ds,label,color,marker in [("QTD000598","Perez NK",GREEN,"o"),("QTD000356","GTEx blood",BLUE,"s"),("QTD000683","Nathan Th1",RED,"D")]:
        rr=sorted([r for r in bfbf if r["dataset_id"]==ds],key=lambda r:float(r["p12"]))
        ins.plot([float(r["p12"]) for r in rr],[float(r["PP.H4"]) for r in rr],color=color,marker=marker,ms=3.5,lw=1.3,label=label)
        inset_endpoints.append((label, color, float(rr[-1]["PP.H4"])))
    ins.axhline(.8,color="#7B8794",ls="--",lw=.7); ins.set_xscale("log"); ins.set_ylim(0,1.03)
    ins.set_xlim(8e-7, 2.8e-5)
    ins.set_xticks([1e-6,1e-5],["1e-6","1e-5"]); ins.tick_params(labelsize=5.7,length=0)
    ins.set_xlabel("p12",fontsize=5.8); ins.set_ylabel("PP.H4",fontsize=5.8)
    label_y = {"Perez NK": .91, "GTEx blood": 1.00, "Nathan Th1": .12}
    for label, color, _ in inset_endpoints:
        ins.text(1.18e-5, label_y[label], label, color=color, fontsize=5.5, va="center", ha="left")
    ins.spines["top"].set_visible(False); ins.spines["right"].set_visible(False)
    ax_d.text(.50,-.23,"Regional compatibility only; not formal coloc.susie",transform=ax_d.transAxes,
              ha="center",fontsize=6.1,color=RED,clip_on=False)
    panel_label(ax_d,"D")

    # E | Allele-aligned direction matrix.
    rows=[("Cordell PBC",-.439362,"7.0e-63","protective"),("FinnGen PBC",-.422314,"2.55e-8","protective"),
          ("UKB PBC",-.408,"6.13e-4","protective"),("OneK1K NK RNA",-.760374,"2.53e-76","lower RNA"),
          ("GTEx blood RNA",-.274502,"5.78e-28","lower RNA"),("CSF protein",.296333,"3.30e-23","higher protein")]
    values=np.array([[r[1]] for r in rows]); dcmap=LinearSegmentedColormap.from_list("dir",[BLUE,"#F7F7F7",RED])
    ax_e.imshow(values,cmap=dcmap,norm=TwoSlopeNorm(vmin=-.8,vcenter=0,vmax=.4),aspect="auto",origin="upper")
    ax_e.set_xlim(-.5,2.6); ax_e.set_ylim(5.5,-.5); ax_e.set_yticks(range(6),[r[0] for r in rows]); ax_e.set_xticks([0,1,2],["Aligned β","P value","Meaning"])
    for i,(_,beta,pv,meaning) in enumerate(rows):
        ax_e.text(0,i,f"{beta:+.3f}",ha="center",va="center",fontsize=7.2,fontweight="bold",color="white" if abs(beta)>.35 else INK)
        ax_e.text(1,i,pv,ha="center",va="center",fontsize=6.3,color=INK)
        ax_e.text(2,i,meaning,ha="center",va="center",fontsize=6.4,fontweight="bold",color=GREEN if beta<0 else RED)
        if i in (2,4): ax_e.axhline(i+.5,color="white",lw=2.5)
    ax_e.set_title("Disease–RNA–protein direction",loc="left",color=INK,pad=8)
    for spine in ax_e.spines.values(): spine.set_visible(False)
    ax_e.tick_params(length=0); ax_e.text(0,-.19,"CSF protein is discordant and compartment-specific.",transform=ax_e.transAxes,fontsize=6.4,color=MUTED)
    panel_label(ax_e,"E")

    fig.suptitle("Figure 3 | Discovery colocalization and external validation of IL12RB2",x=.055,y=.992,ha="left",fontsize=15.2,fontweight="bold",color=INK)
    fig.text(.995,.008,"Discovery sharing, independent molecular replication, FinnGen replication and directional boundaries are shown together.",ha="right",fontsize=6.8,color=MUTED)
    fig.savefig(OUT/"Fig03_IL12RB2_multisignal_final.pdf",pad_inches=.06)
    fig.savefig(OUT/"Fig03_IL12RB2_multisignal_final.png",dpi=300,pad_inches=.06)
    fig.savefig(OUT/"Fig03_IL12RB2_multisignal_final.svg",pad_inches=.06)
    plt.close(fig)


if __name__ == "__main__":
    main()
