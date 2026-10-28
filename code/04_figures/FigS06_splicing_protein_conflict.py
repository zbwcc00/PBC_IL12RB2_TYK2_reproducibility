#!/usr/bin/env python3
"""Supplementary Figure S6: splicing and cross-platform protein conflict."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

REPO=Path(__file__).resolve().parents[2]; ANALYSIS=REPO/"data/derived/analysis"; OUT=REPO/"results/figures/supplementary"; SRC=REPO/"data/derived/figures/supplementary"; OUT.mkdir(parents=True,exist_ok=True); SRC.mkdir(parents=True,exist_ok=True)
AUD=ANALYSIS/"molecular_qtl/il12rb2"
DIR=ANALYSIS/"genetics/causal_direction"
SQ=pd.read_csv(AUD/"16_il12rb2_sqtl_neighboring_splice_audit.tsv",sep="\t"); ISO=pd.read_csv(AUD/"13_il12rb2_transcript_protein_structure.tsv",sep="\t"); MESA=pd.read_csv(AUD/"28_MESA_IL12RB2_shared_cis_variants_aligned.tsv",sep="\t"); COLOC=pd.read_csv(AUD/"31_MESA_IL12RB2_coloc_abf_single_causal_sensitivity.tsv",sep="\t"); CHAIN=pd.read_csv(DIR/"02_allele_direction_chain.tsv",sep="\t"); CONFLICT=pd.read_csv(DIR/"04_protein_splicing_conflict_matrix.tsv",sep="\t")
NAVY="#17324D"; BLUE="#0072B2"; SKY="#56B4E9"; GREEN="#009E73"; AMBER="#E69F00"; VERM="#D55E00"; PURPLE="#8C6BB1"; INK="#26343F"; MID="#687681"; LIGHT="#DCE4E8"; PALE="#F5F8F9"; WHITE="#FFFFFF"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,"axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300,"pdf.fonttype":42})
def header(ax,l,t,s): ax.text(-.03,1.08,l,transform=ax.transAxes,fontsize=13,fontweight="bold",color=WHITE,bbox=dict(boxstyle="round,pad=.23",fc=NAVY,ec="none")); ax.text(.04,1.09,t,transform=ax.transAxes,fontsize=9.4,fontweight="bold",color=NAVY,va="bottom"); ax.text(.04,1.045,s,transform=ax.transAxes,fontsize=6.7,color=MID,va="top")

fig=plt.figure(figsize=(13.2,9.5),facecolor=WHITE); gs=fig.add_gridspec(2,2,left=.075,right=.975,bottom=.09,top=.84,wspace=.28,hspace=.43)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])
fig.suptitle("Supplementary Figure S6 | Splicing and cross-platform protein evidence define a mediation boundary",x=.075,y=.975,ha="left",fontsize=14.2,fontweight="bold",color=NAVY)
fig.text(.075,.937,"Discordant protein assays constrain — rather than invalidate — the RNA-centered locus model.",fontsize=8.5,color=MID)

# A evidence-state matrix from locked conflict table
axA.set_axis_off(); header(axA,"A","Splicing and isoform audit","No validated soluble-IL12RB2 mechanism was identified")
use=CONFLICT.iloc[:3].copy(); state_colors={"unsupported":VERM,"plausible_unvalidated":AMBER,"negative_for_soluble_claim":VERM}
result_display={
 "IL12RB2 sQTL":"No rs6679356–IL12RB2 sQTL in GTEx v10;\ndetected sQTLs belong to SERBP1.",
 "Soluble transcript":"Ensembl transcripts lacking the canonical transmembrane\nsegment are annotation-level candidates only.",
 "Reviewed UniProt isoform":"Reviewed truncated isoforms retain a transmembrane\nsegment and are not validated soluble receptors.",
}
for i,r in use.reset_index(drop=True).iterrows():
 yy=.77-i*.255; c=state_colors.get(r.evidence_state,MID); axA.add_patch(FancyBboxPatch((.03,yy-.13),.94,.20,boxstyle="round,pad=.012,rounding_size=.018",transform=axA.transAxes,fc=PALE,ec=LIGHT)); axA.add_patch(FancyBboxPatch((.055,yy-.065),.24,.075,boxstyle="round,pad=.01",transform=axA.transAxes,fc=c,ec="none")); axA.text(.175,yy-.028,r.question,transform=axA.transAxes,ha="center",va="center",fontsize=7,color=WHITE,fontweight="bold"); axA.text(.33,yy+.015,result_display.get(r.question,r.result),transform=axA.transAxes,fontsize=6.35,color=INK,va="top",linespacing=1.25,clip_on=True); axA.text(.33,yy-.092,r.consequence,transform=axA.transAxes,fontsize=6.25,color=c,fontweight="bold",va="top",clip_on=True)

# B direction chain on harmonized C risk allele
b=CHAIN.copy(); b["beta"]=pd.to_numeric(b.beta_on_risk_allele); labs=["Cordell PBC","FinnGen PBC","UKB PBC","OneK1K NK RNA","GTEx blood RNA","CSF protein"]
cols=[VERM if x>0 else BLUE for x in b.beta]; y=np.arange(len(b)); axB.axvline(0,color=MID,lw=.8); axB.hlines(y,0,b.beta,color=LIGHT,lw=1.2); axB.scatter(b.beta,y,c=cols,s=44,edgecolor=WHITE,lw=.6,zorder=3)
axB.set_yticks(y); axB.set_yticklabels(labs); axB.invert_yaxis(); axB.set_xlim(-.42,.86); axB.set_xlabel("Effect on rs6679356-C (harmonized risk allele)"); axB.grid(axis="x",color=LIGHT,lw=.5)
for i,r in b.reset_index(drop=True).iterrows(): axB.text(r.beta+(.035 if r.beta>=0 else -.035),i,f"{r.beta:+.3f}",ha="left" if r.beta>=0 else "right",va="center",fontsize=6.7,color=INK)
header(axB,"B","Disease–RNA–protein direction chain","CSF protein is opposite to disease and RNA directions")

# C same-sample cross-platform scatter; deterministic density subsampling for raster clarity
x=pd.to_numeric(MESA.z_o,errors="coerce"); yv=pd.to_numeric(MESA.z_s,errors="coerce"); ok=np.isfinite(x)&np.isfinite(yv); x=x[ok].to_numpy(); yv=yv[ok].to_numpy()
lim=np.nanpercentile(np.abs(np.r_[x,yv]),99.5); axC.hexbin(x,yv,gridsize=60,mincnt=1,cmap="Blues",bins="log",extent=(-lim,lim,-lim,lim)); axC.axhline(0,color=MID,lw=.7); axC.axvline(0,color=MID,lw=.7); axC.plot([-lim,lim],[-lim,lim],ls="--",lw=.8,color=VERM); axC.set_xlim(-lim,lim); axC.set_ylim(-lim,lim); axC.set_xlabel("MESA Olink cis-pQTL z"); axC.set_ylabel("MESA SomaScan cis-pQTL z")
pear=np.corrcoef(x,yv)[0,1]; spear=pd.Series(x).corr(pd.Series(yv),method="spearman"); axC.text(.04,.94,f"n={len(x):,}\nPearson r={pear:.3f}\nSpearman ρ={spear:.3f}",transform=axC.transAxes,va="top",fontsize=7.2,color=INK,bbox=dict(boxstyle="round,pad=.35",fc=WHITE,ec=LIGHT))
header(axC,"C","Same-sample cross-platform comparison","8,647 aligned cis variants; density shown on a log scale")

# D coloc sensitivity and locked boundary
cc=COLOC.copy(); cc["label"]=[f"p1=p2={r.p1:g}; p12={r.p12:g}" for _,r in cc.iterrows()]; yy=np.arange(len(cc)); axD.barh(yy,cc["PP.H4"],color=[AMBER,VERM,PURPLE,BLUE],height=.62); axD.set_yticks(yy); axD.set_yticklabels(cc.label); axD.invert_yaxis(); axD.set_xlabel("Olink–SomaScan single-causal PP.H4"); axD.set_xlim(0,.02); axD.grid(axis="x",color=LIGHT,lw=.5)
for i,v in enumerate(cc["PP.H4"]): axD.text(v+.00025,i,f"{v:.4f}",va="center",fontsize=6.8,color=INK)
axD.text(.55,.42,"LOCKED INTERPRETATION\nNo shared primary circulating-protein signal.\nPossible sources: epitope/fragment, PAV or PTM,\nand compartment-specific measurement.",transform=axD.transAxes,fontsize=6.8,color=NAVY,fontweight="bold",bbox=dict(boxstyle="round,pad=.4",fc=PALE,ec=LIGHT))
header(axD,"D","Protein signal-level sensitivity","All tested priors remain non-supportive")

SQ.to_csv(SRC/"FigS06A_sQTL_audit.tsv",sep="\t",index=False); ISO.to_csv(SRC/"FigS06A_isoform_structure.tsv",sep="\t",index=False); b.to_csv(SRC/"FigS06B_allele_direction_chain.tsv",sep="\t",index=False); MESA.to_csv(SRC/"FigS06C_MESA_aligned_cis_variants.tsv",sep="\t",index=False); cc.to_csv(SRC/"FigS06D_crossplatform_coloc.tsv",sep="\t",index=False); CONFLICT.to_csv(SRC/"FigS06_locked_conflict_matrix.tsv",sep="\t",index=False)
fig.text(.075,.025,"Interpretation boundary: CSF or circulating SomaScan/Olink effects cannot be used to infer cell-surface IL12RB2 abundance or pharmacologic direction. These results reject a simple RNA→circulating protein→PBC mediation chain.",fontsize=7,color=VERM)
fig.savefig(OUT/"FigS06_splicing_crossplatform_protein_conflict_final.pdf",bbox_inches="tight",facecolor=WHITE); fig.savefig(OUT/"FigS06_splicing_crossplatform_protein_conflict_final.png",dpi=300,bbox_inches="tight",facecolor=WHITE); print("FigS06 rendered")
