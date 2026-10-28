#!/usr/bin/env python3
"""Supplementary Figure S3: competing genes and locus downgrading.

All quantitative panels are rebuilt from locked project tables. The figure is
descriptive: combinations sharing a PBC locus are not treated as independent.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

REPO = Path(__file__).resolve().parents[2]
ANALYSIS = REPO / "data" / "derived" / "analysis"
OUT = REPO / "results" / "figures" / "supplementary"
SRC = REPO / "data" / "derived" / "figures" / "supplementary"
SRC.mkdir(exist_ok=True)
DISCOVERY = ANALYSIS / "genetics" / "discovery"
REPLICATION = ANALYSIS / "genetics" / "replication"
DISC = pd.read_csv(DISCOVERY / "23a_coloc_susie_combination_summary.tsv", sep="\t")
EXT = pd.read_csv(REPLICATION / "03_external_signal_pair_primary.tsv", sep="\t")
FG = pd.read_csv(REPLICATION / "05_finngen_prespecified_variant_replication.tsv", sep="\t")
REG = pd.read_csv(REPLICATION / "04_finngen_region_qc_top_signal.tsv", sep="\t")

NAVY="#17324D"; BLUE="#0072B2"; SKY="#56B4E9"; GREEN="#009E73"; AMBER="#E69F00"
VERM="#D55E00"; PURPLE="#8C6BB1"; INK="#26343F"; MID="#687681"; LIGHT="#DCE4E8"; PALE="#F5F8F9"; WHITE="#FFFFFF"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,
                     "axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
                     "axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300,"pdf.fonttype":42})

def header(ax, letter, title, sub):
    ax.text(-0.03,1.08,letter,transform=ax.transAxes,fontsize=13,fontweight="bold",color=WHITE,
            bbox=dict(boxstyle="round,pad=.23",fc=NAVY,ec="none"))
    ax.text(0.04,1.09,title,transform=ax.transAxes,fontsize=9.4,fontweight="bold",color=NAVY,va="bottom")
    ax.text(0.04,1.045,sub,transform=ax.transAxes,fontsize=6.7,color=MID,va="top")

fig=plt.figure(figsize=(13.2,9.3),facecolor=WHITE)
gs=fig.add_gridspec(2,2,left=.075,right=.975,bottom=.09,top=.84,wspace=.36,hspace=.43)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])
fig.suptitle("Supplementary Figure S3 | Competing genes and locus-level downgrading",x=.075,y=.975,ha="left",fontsize=15,fontweight="bold",color=NAVY)
fig.text(.075,.937,"Strict-prior signal sharing is separated from transcript resolution and independent disease replication.",fontsize=8.5,color=MID)

# A: locked strict-prior coloc estimates, all requested competing combinations
loci=["2p25.1","1q23.1","18q22.2","17q12"]
genes={"2p25.1":["LINC00299","ID2"],"1q23.1":["FCRL3","FCRL3-AS1"],"18q22.2":["CD226","DOK6"],"17q12":["ORMDL3","GSDMB"]}
colors={"2p25.1":GREEN,"1q23.1":BLUE,"18q22.2":AMBER,"17q12":PURPLE}
rows=[]
for loc in loci:
    for gene in genes[loc]:
        x=DISC[(DISC.pbc_locus==loc)&(DISC.gene_symbol==gene)]
        for _,r in x.iterrows(): rows.append(dict(locus=loc,gene=gene,cell=r.sample_group,pp=float(r["PP.H4"])))
a=pd.DataFrame(rows)
labels=[]; vals=[]; cs=[]
for loc in loci:
    for gene in genes[loc]:
        x=a[(a.locus==loc)&(a.gene==gene)]
        for _,r in x.iterrows():
            labels.append(f"{loc}  {gene} · {r.cell}"); vals.append(r.pp); cs.append(colors[loc])
y=np.arange(len(vals))
axA.hlines(y,0,vals,color=LIGHT,lw=1.2); axA.scatter(vals,y,c=cs,s=42,edgecolor=WHITE,lw=.6,zorder=3)
axA.axvline(.8,color=MID,ls="--",lw=.8); axA.axvline(.9,color=MID,ls=":",lw=.8)
axA.set_yticks(y); axA.set_yticklabels(labels); axA.invert_yaxis(); axA.set_xlim(0,1.02); axA.set_xlabel("Strict-prior coloc.susie PP.H4")
axA.grid(axis="x",color=LIGHT,lw=.5); header(axA,"A","Within-locus competitors","Multiple cell combinations at one locus are correlated evidence")

# B: independent external eQTL support matrix for the two externally audited ambiguous loci
sel=EXT[EXT.locus.isin(["1q23.1","18q22.2"])].copy()
def compact_context(study, group):
    study_map={"Perez_2022":"Perez","Nathan_2022":"Nathan","Randolph_2021":"Randolph",
               "Cytoimmgen":"CytoImmGen","GTEx_v10":"GTEx v10"}
    group_map={"combined_CD4_Naive_UNS_16H":"CD4 naive, 16 h",
               "combined_CD4_Memory_UNS_16H":"CD4 memory, 16 h",
               "CD8+_GZMK+":"CD8 GZMK+","CD4+_Th1":"CD4 Th1","CD4_T_NI":"CD4 T"}
    return study_map.get(str(study),str(study))+" · "+group_map.get(str(group),str(group))
sel["context"]=[compact_context(s,g) for s,g in zip(sel.study_label,sel.sample_group)]
contexts=list(dict.fromkeys(sel.context.tolist()))
glist=["FCRL3","FCRL3-AS1","CD226","DOK6"]
mat=np.full((len(contexts),len(glist)),np.nan)
for i,c in enumerate(contexts):
    for j,g in enumerate(glist):
        z=sel[(sel.context==c)&(sel.gene_symbol==g)]
        if len(z): mat[i,j]=z["PP.H4"].max()
im=axB.imshow(mat,aspect="auto",vmin=0,vmax=1,cmap="viridis")
axB.set_xticks(range(len(glist))); axB.set_xticklabels(glist,rotation=25,ha="right",fontstyle="italic")
axB.set_yticks(range(len(contexts))); axB.set_yticklabels(contexts,fontsize=6.5)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if np.isfinite(mat[i,j]): axB.text(j,i,f"{mat[i,j]:.2f}",ha="center",va="center",fontsize=6.2,color=WHITE if mat[i,j]<.35 or mat[i,j]>.72 else INK)
        else: axB.text(j,i,"—",ha="center",va="center",fontsize=7,color=MID)
cb=fig.colorbar(im,ax=axB,fraction=.035,pad=.025); cb.set_label("PP.H4",fontsize=7)
header(axB,"B","Study-level external eQTL audit","Blank cells mean the gene–context pair was not available")

# C: independent FinnGen layers for audited loci
fg_rows=[]
for loc in ["1p31.3","1q23.1","18q22.2"]:
    q=FG[(FG.locus==loc)&(FG.is_primary_lead.astype(str).str.lower().eq("true"))]
    if len(q) and pd.notna(q.iloc[0].finngen_pvalue_two_sided):
        lead=-np.log10(float(q.iloc[0].finngen_pvalue_two_sided)); lead_lab=f"{q.iloc[0].rsid}"
    else: lead=np.nan; lead_lab="lead unavailable"
    rr=REG[REG.locus==loc].iloc[0]
    fg_rows.append((loc,lead,-np.log10(float(rr.top_pvalue)),-np.log10(float(rr.region_bonferroni_alpha)),lead_lab,bool(rr.region_bonferroni_significant)))
fc=pd.DataFrame(fg_rows,columns=["locus","lead","region","region_threshold","lead_label","region_sig"])
xx=np.arange(len(fc)); w=.28
axC.bar(xx-w/2,fc.lead,w,color=BLUE,label="Prespecified lead"); axC.bar(xx+w/2,fc.region,w,color=AMBER,label="Regional top")
for i,r in fc.iterrows():
    if not np.isfinite(r.lead): axC.text(i-w/2,.25,"QC removed",rotation=90,ha="center",va="bottom",fontsize=6.5,color=VERM)
axC.axhline(-np.log10(.05/3),color=MID,ls="--",lw=.9,label="Lead Bonferroni (.05/3)")
axC.axhline(fc.region_threshold.mean(),color=VERM,ls=":",lw=1.0,label="Regional Bonferroni (~variant-wise)")
axC.set_xticks(xx); axC.set_xticklabels(["IL12RB2\n1p31.3","FCRL3 unit\n1q23.1","CD226/DOK6\n18q22.2"])
axC.set_ylabel("−log10(P), FinnGen"); axC.legend(frameon=False,fontsize=6.6,ncol=2,loc="upper right"); axC.grid(axis="y",color=LIGHT,lw=.5)
header(axC,"C","Independent disease replication","Only IL12RB2 passed prespecified lead and regional replication")

# D: explicit locked decisions
axD.set_axis_off(); header(axD,"D","Evidence lock and downgrade rationale","Gene nomination is bounded by resolution, replication and druggability")
items=[
 ("2p25.1", "LINC00299 > ID2", "Strong contrast", "PP.H4 0.959 vs 0.0004", GREEN, "No validated therapeutic target"),
 ("1q23.1", "FCRL3 / FCRL3-AS1", "Regulatory unit", "External molecular support; FinnGen P=0.321", BLUE, "Sense/antisense transcript unresolved"),
 ("18q22.2", "CD226 / DOK6", "Candidate unit", "PP.H4 0.788–0.797; proxies P=0.296–0.467", AMBER, "Lead QC-removed; causal gene unresolved"),
 ("17q12", "ORMDL3 / GSDMB", "Ambiguous", "PP.H4 0.829 vs 0.738", PURPLE, "Same signal cannot rank causal gene"),
]
for i,(loc,target,grade,evidence,c,reason) in enumerate(items):
    yy=.82-i*.205
    axD.add_patch(FancyBboxPatch((.015,yy-.12),.97,.16,boxstyle="round,pad=.012,rounding_size=.015",transform=axD.transAxes,fc=PALE,ec=LIGHT,lw=.8))
    axD.add_patch(FancyBboxPatch((.03,yy-.082),.13,.085,boxstyle="round,pad=.012",transform=axD.transAxes,fc=c,ec="none"))
    axD.text(.095,yy-.04,loc,transform=axD.transAxes,ha="center",va="center",color=WHITE,fontweight="bold",fontsize=7.4)
    axD.text(.19,yy,target,transform=axD.transAxes,fontweight="bold",fontsize=8.2,color=NAVY)
    axD.text(.19,yy-.047,f"{grade}  |  {evidence}",transform=axD.transAxes,fontsize=6.7,color=INK)
    axD.text(.19,yy-.09,f"Boundary: {reason}",transform=axD.transAxes,fontsize=6.7,color=VERM)

a.to_csv(SRC/"FigS03A_competing_gene_coloc.tsv",sep="\t",index=False)
sel.to_csv(SRC/"FigS03B_external_eqtl_audit.tsv",sep="\t",index=False)
fc.to_csv(SRC/"FigS03C_finngen_replication_layers.tsv",sep="\t",index=False)
pd.DataFrame(items,columns=["locus","target","grade","evidence","color","downgrade_reason"]).drop(columns="color").to_csv(SRC/"FigS03D_evidence_lock.tsv",sep="\t",index=False)
fig.text(.075,.025,"Interpretation boundary: repeated gene–cell combinations at the same locus are not independent discoveries; PP.H4 supports signal sharing, not therapeutic causality.",fontsize=7,color=MID)
fig.savefig(OUT/"FigS03_competing_genes_downgrading_final.pdf",bbox_inches="tight",facecolor=WHITE)
fig.savefig(OUT/"FigS03_competing_genes_downgrading_final.png",dpi=300,bbox_inches="tight",facecolor=WHITE)
print("FigS03 rendered")
