from pathlib import Path
import csv
import gzip
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "figures" / "main"
ANALYSIS = REPO / "data" / "derived" / "analysis"
MESA = ANALYSIS / "molecular_qtl" / "il12rb2"
TYK2 = ANALYSIS / "pharmacology" / "tyk2_safety"
ALIGNED = MESA / "28_MESA_IL12RB2_shared_cis_variants_aligned.tsv"
COLOC = MESA / "31_MESA_IL12RB2_coloc_abf_single_causal_sensitivity.tsv"
META = TYK2 / "genetic" / "02_TYK2_rs34536443_PBC_meta.tsv"
PHEWAS = TYK2 / "genetic" / "01_TYK2_rs34536443_cross_cohort_phewas.tsv"
COND = TYK2 / "locus" / "04_TYK2_two_variant_conditional.tsv"
BOOT = TYK2 / "locus" / "08_TYK2_LD_bootstrap_draws.tsv.gz"

BLUE, SKY, GREEN = "#0072B2", "#56B4E9", "#009E73"
ORANGE, RED, GRAY = "#E69F00", "#D55E00", "#A7B1B8"
INK, MUTED = "#263746", "#667784"


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
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


def main():
    plt.rcParams.update({
        "font.family":"DejaVu Sans", "font.size":8.4, "axes.titlesize":10.2,
        "axes.titleweight":"bold", "axes.labelsize":9, "xtick.labelsize":7.8,
        "ytick.labelsize":7.8, "legend.fontsize":7, "legend.frameon":False,
        "figure.dpi":300, "savefig.dpi":300, "savefig.bbox":"tight",
        "pdf.fonttype":42, "ps.fonttype":42,
    })
    fig=plt.figure(figsize=(13.8,9.1),facecolor="white")
    outer=fig.add_gridspec(2,1,height_ratios=[.92,1.08],hspace=.42)
    top=outer[0].subgridspec(1,2,width_ratios=[1,1.12],wspace=.35)
    bottom=outer[1].subgridspec(1,3,width_ratios=[1.03,1.05,.92],wspace=.40)
    ax_a,ax_b=fig.add_subplot(top[0,0]),fig.add_subplot(top[0,1])
    ax_c,ax_d,ax_e=(fig.add_subplot(bottom[0,i]) for i in range(3))

    aligned=read_tsv(ALIGNED)
    by_key={r["variant_key"]:r for r in aligned}
    selected=[("Olink lead","1:66351133"),("SomaScan lead","1:66560760"),("PBC lead rs6679356","1:67354511")]

    # A | Same-sample cross-platform forest at prespecified/lead variants.
    ybase=np.arange(len(selected))[::-1]*2
    for i,(label,key) in enumerate(selected):
        r=by_key[key]
        for offset,field,sefield,color,marker,plat in [(.22,"beta_o","se_o",BLUE,"o","Olink"),(-.22,"beta_s_aligned","se_s",ORANGE,"s","SomaScan")]:
            beta,se=float(r[field]),float(r[sefield]); y=ybase[i]+offset
            ax_a.errorbar(beta,y,xerr=1.96*se,fmt=marker,ms=5.5,capsize=2.5,color=color,lw=1.3,
                          label=plat if i==0 else None)
    ax_a.axvline(0,color="#6B737A",lw=.9)
    ax_a.set_yticks(ybase,[x[0] for x in selected]); ax_a.set_ylim(-1,5)
    ax_a.set_xlabel("Allele-aligned cis-pQTL effect (assay units)")
    ax_a.set_title("Same-sample Olink–SomaScan pQTL comparison",loc="left",color=INK,pad=8)
    ax_a.legend(loc="lower left"); style_axis(ax_a,"x")
    ax_a.text(.01,-.19,"MESA n=1,889; effect scales are platform-specific.",transform=ax_a.transAxes,fontsize=6.8,color=MUTED)
    panel_label(ax_a,"A")

    # B | Genome-wide cis-window cross-platform z-score comparison.
    z_o=np.array([float(r["z_o"]) for r in aligned]); z_s=np.array([float(r["z_s"]) for r in aligned])
    minp=np.array([float(r["min_p"]) for r in aligned]); sig=minp<.01
    ax_b.scatter(z_o[~sig],z_s[~sig],s=5,color="#CAD2D7",alpha=.35,linewidths=0,rasterized=True)
    ax_b.scatter(z_o[sig],z_s[sig],s=10,color=SKY,alpha=.7,linewidths=0,rasterized=True,label="Either-platform P<0.01")
    for label,key,color in [("Olink lead","1:66351133",BLUE),("SomaScan lead","1:66560760",ORANGE),("rs6679356","1:67354511",GREEN)]:
        r=by_key[key]; xo,ys=float(r["z_o"]),float(r["z_s"])
        ax_b.scatter(xo,ys,s=50,color=color,edgecolor="white",lw=.7,zorder=4)
        ax_b.annotate(label,(xo,ys),xytext=(5,6),textcoords="offset points",fontsize=6.6,color=color)
    lim=max(4,np.nanpercentile(np.abs(np.r_[z_o,z_s]),99.5)); ax_b.plot([-lim,lim],[-lim,lim],ls="--",lw=.8,color="#8A939A")
    ax_b.axhline(0,color="#AAB2B8",lw=.6); ax_b.axvline(0,color="#AAB2B8",lw=.6)
    ax_b.set_xlim(-lim,lim); ax_b.set_ylim(-lim,lim)
    ax_b.set_xlabel("Olink cis-pQTL z score"); ax_b.set_ylabel("SomaScan cis-pQTL z score")
    ax_b.set_title("Cross-platform cis-window concordance is weak",loc="left",color=INK,pad=8)
    ax_b.legend(loc="lower right"); style_axis(ax_b)
    coloc=read_tsv(COLOC); max_h4=max(float(r["PP.H4"]) for r in coloc)
    ax_b.text(.02,.97,"8,647 aligned variants\nDirection agreement 47.45%\nPearson r = −0.084\nABF PP.H4 ≤ 0.0166",
              transform=ax_b.transAxes,ha="left",va="top",fontsize=7,color=RED,
              bbox=dict(facecolor="white",edgecolor="#E2E6E9",boxstyle="round,pad=.3",alpha=.92))
    panel_label(ax_b,"B")

    # C | TYK2 partial loss-of-function meta-analysis.
    phe=[r for r in read_tsv(PHEWAS) if r["domain"]=="PBC"]
    meta=read_tsv(META)[0]
    labels=["FinnGen R12","UKB–TOPMed","Fixed-effect meta"]
    effects=[float(r["odds_ratio"]) for r in phe]+[float(meta["odds_ratio"])]
    lows=[]; highs=[]; pvals=[]
    for r in phe:
        b,se=float(r["beta_per_C_allele"]),float(r["se"]); lows.append(math.exp(b-1.96*se)); highs.append(math.exp(b+1.96*se)); pvals.append(float(r["p_value"]))
    lows.append(float(meta["ci95_low"])); highs.append(float(meta["ci95_high"])); pvals.append(float(meta["p_value"]))
    yy=np.arange(3)[::-1]
    for y,label,e,lo,hi,pv,color,marker in zip(yy,labels,effects,lows,highs,pvals,[BLUE,SKY,GREEN],["o","s","D"]):
        ax_c.errorbar(e,y,xerr=[[e-lo],[hi-e]],fmt=marker,ms=6,capsize=3,color=color,lw=1.4)
        ax_c.text(1.02,y,f"OR {e:.2f} [{lo:.2f}, {hi:.2f}]\nP={pv:.2e}",va="center",fontsize=6.5,color=color)
    ax_c.axvline(1,color="#6B737A",lw=.9); ax_c.set_xlim(.25,1.52); ax_c.set_ylim(-.55,2.55)
    ax_c.set_yticks(yy,labels); ax_c.set_xlabel("PBC odds ratio per rs34536443-C allele")
    ax_c.set_title("TYK2 partial loss-of-function meta-analysis",loc="left",color=INK,pad=8)
    style_axis(ax_c,"x")
    ax_c.text(.02,.04,"Fixed effect: OR=0.560; I²=0%\n1,016 cases / 770,550 controls",transform=ax_c.transAxes,fontsize=6.8,color=MUTED)
    panel_label(ax_c,"C")

    # D | Two-variant conditional approximation.
    cond=[r for r in read_tsv(COND) if r["dataset"]=="GCST005581" and r["LD_definition"]=="unphased_dosage_r_primary"]
    names=["rs34536443-C","rs2304256-A"]
    marginal=[float(r["marginal_z"]) for r in cond]; conditional=[float(r["conditional_z_approx"]) for r in cond]
    y=np.arange(2)[::-1]
    for yy,name,m,c in zip(y,names,marginal,conditional):
        ax_d.plot([m,c],[yy,yy],color="#B9C1C7",lw=1.6)
        ax_d.scatter(m,yy,s=52,color=GRAY,marker="o",label="Marginal" if yy==1 else None,zorder=3)
        ax_d.scatter(c,yy,s=58,color=BLUE,marker="D",label="Conditional" if yy==1 else None,zorder=4)
        p=2*norm.sf(abs(c)); ax_d.text(c+.16,yy-.16,f"P={p:.2e}",fontsize=6.6,color=BLUE)
    ax_d.axvline(0,color="#6B737A",lw=.9); ax_d.set_xlim(-7.8,.4); ax_d.set_ylim(-.6,1.6)
    ax_d.set_yticks(y,names); ax_d.set_xlabel("PBC z score")
    ax_d.set_title("Two-variant conditional approximation",loc="left",color=INK,pad=8)
    ax_d.legend(loc="upper right"); style_axis(ax_d,"x")
    ax_d.text(.02,.04,"1000G EUR signed dosage r=0.328\nSummary-z approximation; not individual-level COJO",transform=ax_d.transAxes,fontsize=6.6,color=RED)
    panel_label(ax_d,"D")

    # E | 100,000-draw LD bootstrap distributions.
    draws=read_tsv(BOOT)
    scenarios=[("EUR_all_population_stratified","EUR all",BLUE),("FIN","FIN",ORANGE),("EUR_non_FIN_population_stratified","EUR non-FIN",GREEN)]
    data=[]; labels=[]; colors=[]
    for sc,label,color in scenarios:
        vals=[float(r["signed_r"]) for r in draws if r["scenario"]==sc and r["valid"]=="True"]
        data.append(vals); labels.append(label); colors.append(color)
    vp=ax_e.violinplot(data,positions=np.arange(3),showmedians=True,showextrema=False,widths=.78)
    for body,color in zip(vp["bodies"],colors): body.set_facecolor(color); body.set_edgecolor("white"); body.set_alpha(.72)
    vp["cmedians"].set_color(INK); vp["cmedians"].set_linewidth(1.3)
    ax_e.set_xticks(range(3),labels,rotation=12); ax_e.set_ylabel("Signed LD r (effect alleles)")
    ax_e.set_title("100,000-draw LD bootstrap",loc="left",color=INK,pad=8)
    style_axis(ax_e,"y")
    ax_e.text(.03,.97,"Both conditional P<0.05:\nEUR all 100%\nFIN 99.998%\nEUR non-FIN 100%",
              transform=ax_e.transAxes,ha="left",va="top",fontsize=6.7,color=INK,
              bbox=dict(facecolor="white",edgecolor="#E2E6E9",boxstyle="round,pad=.25",alpha=.92))
    ax_e.text(.02,.03,"FIN rare-allele valid draws: 95,233/100,000",transform=ax_e.transAxes,fontsize=6.3,color=MUTED)
    panel_label(ax_e,"E")

    fig.suptitle("Figure 4 | Cross-platform protein boundary and TYK2 human genetic evidence",x=.055,y=.992,ha="left",fontsize=15.2,fontweight="bold",color=INK)
    fig.text(.995,.008,"Protein assays do not replicate a shared IL12RB2 signal; TYK2 evidence supports a distinct, genetically anchored drug node.",ha="right",fontsize=6.8,color=MUTED)
    fig.savefig(OUT/"Fig04_pQTL_TYK2_final.pdf",pad_inches=.06)
    fig.savefig(OUT/"Fig04_pQTL_TYK2_final.png",dpi=300,pad_inches=.06)
    plt.close(fig)


if __name__=="__main__":
    main()
