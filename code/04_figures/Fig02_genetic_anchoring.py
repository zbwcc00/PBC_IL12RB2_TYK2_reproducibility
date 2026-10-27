from pathlib import Path
import gzip
import math
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "figures" / "main"
SOURCE = REPO / "data" / "derived" / "figures" / "main" / "Fig02"
RAW_GWAS = SOURCE / "Fig02A_manhattan_thinned.tsv.gz"
SIGNIFICANT = SOURCE / "Fig02A_genomewide_significant.tsv.gz"
LEADS = SOURCE / "Fig02A_lead_clusters.tsv"
LOCUS = SOURCE / "Fig02B_IL12RB2_locus_with_lead_r2.tsv.gz"
CS = SOURCE / "Fig02C_credible_sets.tsv"
ALPHA = SOURCE / "Fig02D_susie_alpha.tsv"

INK = "#253746"
MUTED = "#687984"
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
AMBER = "#E69F00"
VERM = "#D55E00"
PINK = "#CC79A7"
GRAY = "#B8C1C8"
LIGHT = "#EDF2F5"
SIGCOL = {1: BLUE, 2: AMBER, 3: GREEN}
CHR_COLORS = [
    "#42A5C8", "#F2C230", "#443A91", "#566136", "#D41465", "#F07A23",
    "#C477A3", "#7FCDBB", "#9B192E", "#1E7F63", "#4396C1", "#F1C52E",
    "#42358B", "#606534", "#D31363", "#F27623", "#C778A4", "#77C8B6",
    "#981B2D", "#26836A", "#4A9EC4", "#F4C938",
]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.0,
    "axes.titlesize": 9.3, "axes.titleweight": "bold",
    "axes.labelsize": 8.2, "xtick.labelsize": 7.0, "ytick.labelsize": 7.0,
    "legend.fontsize": 6.8, "legend.frameon": False,
    "figure.dpi": 300, "savefig.dpi": 300,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def panel(ax, letter, title):
    ax.text(-.10, 1.06, letter, transform=ax.transAxes, fontsize=12.5,
            fontweight="bold", color=INK, va="top")
    ax.set_title(title, loc="left", pad=6, color=INK)


def clean(ax, grid="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis=grid, color="#DDE5EA", lw=.55, alpha=.7)
        ax.set_axisbelow(True)


def load_gwas():
    return pd.read_csv(RAW_GWAS, sep="	")


def chromosome_offsets(d):
    maxpos=d.groupby("chromosome").base_pair_location.max().reindex(range(1,23))
    offsets={}; cur=0
    for chrom in range(1,23):
        offsets[chrom]=cur
        cur += int(maxpos.loc[chrom]) + 12_000_000
    return offsets, maxpos


def read_ld_lead(rsids, lead):
    locus = pd.read_csv(LOCUS, sep="	")
    lookup = dict(zip(locus["rsid"].astype(str), locus["lead_r2"]))
    return np.array([lookup.get(str(rsid), np.nan) for rsid in rsids])


def genes_track(ax, xmin, xmax):
    genes=[
        (67.465,67.520,"SLC35D1",0), (67.558,67.698,"C1orf141",1),
        (67.632,67.726,"IL23R",0), (67.773,67.863,"IL12RB2",1),
        (67.873,67.896,"SERBP1",0), (68.151,68.154,"GADD45A",1),
        (68.167,68.299,"GNG12",0), (68.512,68.517,"DIRAS3",1),
    ]
    ax.set_xlim(xmin,xmax); ax.set_ylim(-.3,2.1)
    ax.set_yticks([])
    for side in ("top","left","right"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#93A2AC")
    ax.tick_params(axis="x", labelsize=6.5, colors=MUTED, length=2.5)
    for a,b,name,row in genes:
        if b<xmin or a>xmax: continue
        y=.45+row*.80
        a=max(a,xmin); b=min(b,xmax)
        ax.plot([a,b],[y,y],color=INK,lw=1.2)
        step=max((b-a)/7,.008)
        for xx in np.arange(a,b,step):
            ax.annotate("",(min(xx+step*.8,b),y),(xx,y),arrowprops=dict(arrowstyle="-|>",lw=.7,color=INK))
        ax.text((a+b)/2,y+.20,name,ha="center",va="bottom",fontsize=6.1,
                color=VERM if name=="IL12RB2" else MUTED,
                fontweight="bold" if name=="IL12RB2" else "normal")


def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    fig=plt.figure(figsize=(12.8,9.0), facecolor="white")
    gs=fig.add_gridspec(3,3,height_ratios=[.95,1.22,.95],width_ratios=[1.16,1.16,.94],
                        left=.070,right=.970,top=.925,bottom=.080,wspace=.36,hspace=.45)
    subA=gs[0,:].subgridspec(2,1,height_ratios=[11,1.15],hspace=.025)
    axA=fig.add_subplot(subA[0]); axAr=fig.add_subplot(subA[1],sharex=axA)
    subB=gs[1,:2].subgridspec(2,1,height_ratios=[4.5,1],hspace=.04)
    axB=fig.add_subplot(subB[0]); axBg=fig.add_subplot(subB[1],sharex=axB)
    axC=fig.add_subplot(gs[1,2])
    axD=fig.add_subplot(gs[2,:])

    # A — genome-wide discovery background
    g=load_gwas(); offsets,maxpos=chromosome_offsets(g)
    g["gx"]=[p+offsets[int(c)] for c,p in zip(g.chromosome,g.base_pair_location)]
    g["mlogp"]=-np.log10(np.clip(g.p_value,1e-300,None))
    for chrom,dd in g.groupby("chromosome"):
        axA.scatter(dd.gx,dd.mlogp,s=2.3,c=CHR_COLORS[int(chrom)-1],alpha=.78,lw=0,rasterized=True)
    # Restricted labels: manuscript anchors and a few independent reference peaks.
    lead=pd.read_csv(LEADS,sep="\t")
    label_map={
        "rs7774434":("HLA", "#6C7A84"),
        "rs6679356":("IL12RB2", BLUE),
        "rs485789":("rs485789", VERM),
        "rs7674640":("NFKB1", VERM),
        "rs33938760":("GSDMB/ORMDL3", VERM),
        "rs2304256":("TYK2", "#D55E00"),
    }
    vertical={"rs7774434":10,"rs6679356":9,"rs485789":7,"rs7674640":8,"rs33938760":8,"rs2304256":8}
    for rsid,(name,color) in label_map.items():
        r=lead[lead.variant_id.eq(rsid)].iloc[0]
        xx=offsets[int(r.chromosome)]+int(r.base_pair_location); yy=-math.log10(float(r.p_value))
        axA.scatter(xx,yy,s=29,c=color,ec="white",lw=.65,zorder=5)
        axA.annotate(f"{name}\n{rsid}",(xx,yy),xytext=(xx,min(yy+vertical[rsid],110)),textcoords="data",
                     ha="center",va="bottom",fontsize=5.9,
                     fontweight="bold" if rsid in ("rs6679356","rs2304256") else "normal",color=color,
                     arrowprops=dict(arrowstyle="-",color=color,lw=.55,shrinkA=1,shrinkB=2))
    axA.axhline(-math.log10(5e-8),color="#7D8992",ls="--",lw=.8)
    centers=[offsets[c]+maxpos.loc[c]/2 for c in range(1,23)]
    axA.set_xticks([])
    axA.set_xlim(-5e6,max(g.gx)+5e6); axA.set_ylim(0,max(112,g.mlogp.max()*1.06))
    axA.set_ylabel(r"$-\log_{10}(P)$")
    panel(axA,"A","PBC genome-wide association landscape")
    clean(axA,None)
    axA.text(.995,.96,"GCST90061440  •  European ancestry  •  GRCh37",transform=axA.transAxes,
             ha="right",va="top",fontsize=6.6,color=MUTED)
    # Significant-variant rug follows the supplied reference style.
    sig=pd.read_csv(SIGNIFICANT,sep="\t")
    sigx=[offsets[int(c)]+int(p) for c,p in zip(sig.chromosome,sig.base_pair_location)]
    axAr.vlines(sigx,0,1,color="#207A3B",lw=.20,alpha=.68)
    axAr.set_ylim(0,1); axAr.set_yticks([]); axAr.set_xticks(centers); axAr.set_xticklabels(range(1,23))
    axAr.set_xlabel("Chromosome",labelpad=2)
    for side in ("top","left","right"): axAr.spines[side].set_visible(False)
    axAr.spines["bottom"].set_color("#84939D")
    axAr.tick_params(axis="x",length=2.5,pad=1)
    axAr.text(-.012,-.08,"Chr",transform=axAr.transAxes,ha="right",va="top",fontsize=6.5,fontweight="bold")
    axAr.text(.995,.66,"Genome-wide significant variants",transform=axAr.transAxes,ha="right",fontsize=5.8,color="#207A3B")

    # B — regional association + LD + gene track
    loc=pd.read_csv(LOCUS,sep="\t")
    loc["mb"]=loc.position_grch37_1based/1e6
    loc["mlogp"]=-np.log10(np.clip(loc.pbc_pvalue,1e-300,None))
    loc["r2"]=read_ld_lead(loc.rsid.astype(str).tolist(),"rs6679356")
    cmap=LinearSegmentedColormap.from_list("ld",["#D8DDE1",SKY,BLUE,VERM])
    sc=axB.scatter(loc.mb,loc.mlogp,c=loc.r2,cmap=cmap,norm=Normalize(0,1),s=11,lw=0,alpha=.88,rasterized=True)
    cs=pd.read_csv(CS,sep="\t")
    cs=cs[(cs.locus_id=="PBC04")&(cs.model_L==10)].copy()
    csmap={}
    for _,r in cs.iterrows():
        for snp in str(r.credible_set_snps).split(";"): csmap[snp]=int(r.susie_component)
    for sig,color in SIGCOL.items():
        d=loc[loc.rsid.map(csmap).eq(sig)]
        axB.scatter(d.mb,d.mlogp,s=40,facecolor=color,edgecolor="white",lw=.7,zorder=4,label=f"CS{sig}")
    lead=loc[loc.rsid.eq("rs6679356")].iloc[0]
    axB.annotate("rs6679356",(lead.mb,lead.mlogp),xytext=(8,10),textcoords="offset points",
                 fontsize=7.2,fontweight="bold",color=BLUE)
    axB.axhline(-math.log10(5e-8),color="#7D8992",ls="--",lw=.8)
    axB.set_ylabel(r"$-\log_{10}(P)$"); axB.tick_params(labelbottom=False)
    axB.legend(loc="upper right",ncol=3,title="95% credible set",title_fontsize=6.6)
    panel(axB,"B","IL12RB2 regional association, lead-variant LD and gene context")
    clean(axB,None)
    cbar=fig.colorbar(sc,ax=axB,fraction=.022,pad=.012); cbar.set_label(r"$r^2$ with rs6679356",fontsize=7); cbar.ax.tick_params(labelsize=6)
    genes_track(axBg,loc.mb.min(),loc.mb.max()); axBg.set_xlabel("Position on chromosome 1 (Mb; GRCh37)",fontsize=8)

    # C — credible-set composition (conditional inclusion probabilities)
    alpha=pd.read_csv(ALPHA,sep="\t")
    comp=[]
    for _,r in cs.sort_values("susie_component").iterrows():
        sig=int(r.susie_component); col=f"signal{sig}"
        snps=str(r.credible_set_snps).split(";")
        q=alpha[alpha.rsid.isin(snps)][["rsid",col]].sort_values(col,ascending=True)
        for _,v in q.iterrows(): comp.append((sig,v.rsid,float(v[col])))
    y=np.arange(len(comp)); vals=[v[2] for v in comp]; labels=[v[1] for v in comp]
    axC.barh(y,vals,color=[SIGCOL[v[0]] for v in comp],height=.66)
    axC.set_yticks(y); axC.set_yticklabels(labels,fontsize=6.4); axC.set_xlim(0,1)
    axC.set_xlabel("Conditional inclusion probability")
    for yy,(sig,snp,val) in zip(y,comp):
        axC.text(min(val+.025,.96),yy,f"{val:.2f}",va="center",fontsize=6.2,color=INK)
    panel(axC,"C","95% credible-set composition")
    clean(axC,"x")
    axC.text(.98,.03,"L=10: CS sizes 2, 3 and 2",transform=axC.transAxes,ha="right",fontsize=6.5,color=MUTED)

    # D — three SuSiE component tracks across the region
    m=loc[["rsid","mb"]].merge(alpha,on="rsid",how="inner").sort_values("mb")
    for sig,offset in zip((1,2,3),(2,1,0)):
        a=m[f"signal{sig}"].to_numpy()
        axD.scatter(m.mb,np.full(len(m),offset)+a*.72,s=3.5,color=SIGCOL[sig],alpha=.18,lw=0,rasterized=True)
        top=m.nlargest(4,f"signal{sig}")
        axD.scatter(top.mb,offset+top[f"signal{sig}"]*.72,s=18,color=SIGCOL[sig],ec="white",lw=.45,zorder=4)
        leadrow=top.iloc[0]
        axD.text(leadrow.mb,offset+float(leadrow[f"signal{sig}"])*.72+.08,str(leadrow.rsid),
                 ha="center",fontsize=5.8,color=SIGCOL[sig])
    axD.set_yticks([0,1,2]); axD.set_yticklabels(["Signal 3","Signal 2","Signal 1"])
    axD.set_ylim(-.18,2.92); axD.set_xlabel("Position on chromosome 1 (Mb)")
    axD.set_ylabel("SuSiE component")
    panel(axD,"D","susie_rss resolves three PBC signals")
    clean(axD,None)
    axD.text(.98,.04,"L=10  •  3,638 variants  •  converged",transform=axD.transAxes,ha="right",fontsize=6.4,color=MUTED)

    fig.suptitle("Figure 2 | PBC GWAS anchoring and fine-mapping of the IL12RB2 locus",
                 x=.065,y=.982,ha="left",fontsize=15,fontweight="bold",color=INK)
    fig.text(.975,.025,
             "TYK2 is marked as a genome-wide disease locus; locus-level fine-mapping panels refer to IL12RB2. LD: 1000 Genomes EUR.",
             ha="right",fontsize=6.6,color=MUTED)
    for stem in ("Fig02_genetic_anchoring_final",):
        fig.savefig(OUT/f"{stem}.pdf",bbox_inches="tight",pad_inches=.06)
        fig.savefig(OUT/f"{stem}.png",dpi=300,bbox_inches="tight",pad_inches=.06)
    plt.close(fig)


if __name__=="__main__":
    main()
