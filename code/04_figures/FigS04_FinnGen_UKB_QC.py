#!/usr/bin/env python3
"""Supplementary Figure S4: multilevel IL12RB2 disease replication and method boundary."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

REPO=Path(__file__).resolve().parents[2]; ANALYSIS=REPO/"data/derived/analysis"; OUT=REPO/"results/figures/supplementary"; SRC=REPO/"data/derived/figures/supplementary"; OUT.mkdir(parents=True,exist_ok=True); SRC.mkdir(parents=True,exist_ok=True)
BASE=ANALYSIS/"genetics/replication"
DIR=ANALYSIS/"genetics/causal_direction"
CHAIN=pd.read_csv(DIR/"02_allele_direction_chain.tsv",sep="\t")
FG=pd.read_csv(BASE/"05_finngen_prespecified_variant_replication.tsv",sep="\t")
REG=pd.read_csv(BASE/"04_finngen_region_qc_top_signal.tsv",sep="\t")
CS=pd.read_csv(BASE/"16_il12rb2_finngen_single_causal_abf_95cs.tsv",sep="\t")
BF=pd.read_csv(BASE/"17_il12rb2_finngen_external_eqtl_bf_bf_sensitivity.tsv",sep="\t")
OV=pd.read_csv(BASE/"19_il12rb2_finngen_external_eqtl_95cs_overlap.tsv",sep="\t")
LOCK=pd.read_csv(ANALYSIS/"evidence/final_evidence_lock.tsv",sep="\t")

NAVY="#17324D"; BLUE="#0072B2"; SKY="#56B4E9"; GREEN="#009E73"; AMBER="#E69F00"; VERM="#D55E00"; PURPLE="#8C6BB1"; INK="#26343F"; MID="#687681"; LIGHT="#DCE4E8"; PALE="#F5F8F9"; WHITE="#FFFFFF"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,"axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300,"pdf.fonttype":42})
def header(ax,l,t,s):
 ax.text(-.03,1.08,l,transform=ax.transAxes,fontsize=13,fontweight="bold",color=WHITE,bbox=dict(boxstyle="round,pad=.23",fc=NAVY,ec="none")); ax.text(.04,1.09,t,transform=ax.transAxes,fontsize=9.4,fontweight="bold",color=NAVY,va="bottom"); ax.text(.04,1.045,s,transform=ax.transAxes,fontsize=6.7,color=MID,va="top")

fig=plt.figure(figsize=(13.2,9.6),facecolor=WHITE); gs=fig.add_gridspec(2,2,left=.075,right=.975,bottom=.09,top=.84,wspace=.30,hspace=.45)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])
fig.suptitle("Supplementary Figure S4 | Independent disease replication and regional compatibility at IL12RB2",x=.075,y=.975,ha="left",fontsize=14.5,fontweight="bold",color=NAVY)
fig.text(.075,.937,"Allele-aligned cohort replication is distinguished from single-causal BF–BF compatibility.",fontsize=8.5,color=MID)

# A risk allele C, SE reconstructed for UKB from p and beta when unavailable
d=CHAIN[CHAIN.layer.str.contains("PBC GWAS",na=False)].copy()
d["beta"]=pd.to_numeric(d.beta_on_risk_allele); d["p"]=pd.to_numeric(d.p_value)
from scipy.stats import norm
def p_to_z(p): return norm.isf(p/2)
d["se"]=[FG.loc[(FG.locus=="1p31.3")&(FG.rsid=="rs6679356"),"finngen_se"].iloc[0] if "FinnGen" in s else abs(b)/p_to_z(p) for s,b,p in zip(d.source,d.beta,d.p)]
d["label"]=["Cordell discovery","FinnGen R12","UK Biobank"]
y=np.arange(len(d)); lo=d.beta-1.96*d.se; hi=d.beta+1.96*d.se
axA.errorbar(d.beta,y,xerr=[d.beta-lo,hi-d.beta],fmt="o",color=BLUE,ecolor=BLUE,capsize=3,lw=1.4)
axA.axvline(0,color=MID,lw=.8); axA.set_yticks(y); axA.set_yticklabels(d.label); axA.invert_yaxis(); axA.set_xlabel("Log-odds effect of rs6679356-C (risk allele)")
for i,r in d.reset_index(drop=True).iterrows(): axA.text(hi.iloc[i]+.018,i,f"β={r.beta:.3f}; P={r.p:.2g}",va="center",fontsize=6.8,color=INK)
axA.set_xlim(0.15,.63); axA.grid(axis="x",color=LIGHT,lw=.5); header(axA,"A","Allele-aligned replication","C increased PBC risk in all three disease datasets")

# B multilevel ladder
axB.set_axis_off(); header(axB,"B","FinnGen replication ladder","Variant, region and gene tests answer different questions")
gene=LOCK[LOCK.evidence_id=="E013"].iloc[0]
levels=[("Prespecified variant","rs6679356",2.54935e-8,"Bonferroni replicated",GREEN),
        ("Regional top","rs6659932",float(REG.loc[REG.locus=="1p31.3","top_pvalue"].iloc[0]),"Genome-wide significant",GREEN),
        ("Gene-level","IL12RB2",1.81e-5,"BH-FDR 0.00484; Bonferroni P=0.440",AMBER)]
for i,(level,item,p,status,c) in enumerate(levels):
 yy=.78-i*.25; axB.add_patch(FancyBboxPatch((.04,yy-.105),.92,.16,boxstyle="round,pad=.012,rounding_size=.018",transform=axB.transAxes,fc=PALE,ec=LIGHT)); axB.add_patch(FancyBboxPatch((.065,yy-.065),.18,.08,boxstyle="round,pad=.012",transform=axB.transAxes,fc=c,ec="none")); axB.text(.155,yy-.025,level,transform=axB.transAxes,ha="center",va="center",color=WHITE,fontweight="bold",fontsize=7); axB.text(.285,yy,item,transform=axB.transAxes,fontweight="bold",fontsize=8.2,color=NAVY); axB.text(.285,yy-.048,f"P={p:.2g}  |  {status}",transform=axB.transAxes,fontsize=7,color=INK)

# C 95% single-causal credible set at prior SD=.2
c=CS[CS.prior_sd_key.astype(str)=="0.2"].copy().sort_values("rank")
top=c.head(12).copy(); rest=max(0,1-top.posterior.sum()); labs=[str(x)[:18] for x in top.rsids]; vals=top.posterior.tolist()
if rest>0: labs.append(f"Other {len(c)-len(top)} CS variants"); vals.append(rest)
y=np.arange(len(vals)); axC.barh(y,vals,color=[BLUE]*min(12,len(vals))+([LIGHT] if len(vals)>12 else []),height=.65)
axC.set_yticks(y); axC.set_yticklabels(labs); axC.invert_yaxis(); axC.set_xlabel("FinnGen regional posterior (Wakefield ABF)"); axC.grid(axis="x",color=LIGHT,lw=.5)
header(axC,"C","FinnGen single-causal 95% credible set",f"{len(c)} variants at prior SD 0.20; regional top rs6659932")

# D sensitivity plus CS overlap
order=["Perez_2022","GTEx_v10","Nathan_2022"]; names={"Perez_2022":"Perez NK","GTEx_v10":"GTEx blood","Nathan_2022":"Nathan Th1"}; pri=[1e-6,5e-6,1e-5]; offs=[-.20,0,.20]; cols=[BLUE,AMBER,GREEN]
for off,pv,col in zip(offs,pri,cols):
 x=[]
 for st in order:
  z=BF[(BF.study_label==st)&np.isclose(BF.p12,pv)]
  x.append(float(z.iloc[0]["PP.H4"]) if len(z) else np.nan)
 axD.scatter(x,np.arange(3)+off,s=38,color=col,label=f"p12={pv:g}",edgecolor=WHITE,lw=.5,zorder=3)
for i,st in enumerate(order):
 ov=OV[OV.study_label==st].iloc[0]; axD.text(1.015,i,f"CS overlap {int(ov.shared_95cs_variant_n)}/{int(ov.external_95cs_variant_n)}",va="center",fontsize=6.5,color=INK)
axD.axvline(.8,color=MID,ls="--",lw=.8); axD.set_yticks(range(3)); axD.set_yticklabels([names[x] for x in order]); axD.invert_yaxis(); axD.set_xlim(0,1.27); axD.set_xlabel("FinnGen–eQTL BF–BF PP.H4"); axD.grid(axis="x",color=LIGHT,lw=.5); axD.legend(frameon=False,fontsize=6.5,loc="upper left")
header(axD,"D","Prior sensitivity and credible-set geometry","Nathan Th1 is a prespecified informative negative context")

d.to_csv(SRC/"FigS04A_allele_aligned_GWAS.tsv",sep="\t",index=False); pd.DataFrame(levels,columns=["level","item","pvalue","status","color"]).drop(columns="color").to_csv(SRC/"FigS04B_replication_ladder.tsv",sep="\t",index=False); c.to_csv(SRC/"FigS04C_FinnGen_ABF_95CS.tsv",sep="\t",index=False); BF.to_csv(SRC/"FigS04D_BF_BF_sensitivity.tsv",sep="\t",index=False); OV.to_csv(SRC/"FigS04D_credible_set_overlap.tsv",sep="\t",index=False)
fig.text(.075,.025,"Method boundary: BF–BF is a single-causal approximation. Complete, unthresholded Finnish signed LD was unavailable; therefore this is not FinnGen coloc.susie. The public 19-variant official matrix was used only for structural checks.",fontsize=7,color=VERM)
fig.savefig(OUT/"FigS04_FinnGen_UKB_replication_QC_final.pdf",bbox_inches="tight",facecolor=WHITE); fig.savefig(OUT/"FigS04_FinnGen_UKB_replication_QC_final.png",dpi=300,bbox_inches="tight",facecolor=WHITE); print("FigS04 rendered")
