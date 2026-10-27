from pathlib import Path as FilePath
import math

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle, PathPatch
from matplotlib.path import Path as MplPath


OUT = FilePath(__file__).resolve().parents[2] / "results" / "figures" / "main"

INK = "#233446"
MUTED = "#647688"
BG = "#F7FAFC"
WHITE = "#FFFFFF"
LINE = "#DCE6EC"
BLUE = "#1683BD"
TEAL = "#12A186"
AMBER = "#ECA538"
VIOLET = "#7969B2"
ORANGE = "#E87336"
GREEN = "#299D68"
RED = "#D9534F"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def rr(ax, x, y, w, h, fc=WHITE, ec=LINE, lw=1, radius=.015, z=2):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0.004,rounding_size={radius}",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z)
    ax.add_patch(p)
    return p


def txt(ax, x, y, s, size=8, color=INK, weight="normal", ha="left", va="top", z=7, spacing=1.22):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
            ha=ha, va=va, zorder=z, linespacing=spacing)


def badge(ax, x, y, w, h, s, fc, color, size=6.7, weight="bold"):
    rr(ax, x, y, w, h, fc=fc, ec="none", radius=h/2, z=5)
    txt(ax, x+w/2, y+h/2, s, size, color, weight, "center", "center", 6)


def arrow(ax, x1, x2, y=.575, color="#91A2AE"):
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                                 mutation_scale=15, color=color, linewidth=1.7,
                                 shrinkA=0, shrinkB=0, zorder=9))


def card(ax, x, y, w, h, color, n, title):
    rr(ax, x+.004, y-.006, w, h, fc="#D7E2E9", ec="none", radius=.020, z=1)
    rr(ax, x, y, w, h, fc=WHITE, ec=LINE, lw=.9, radius=.020, z=2)
    # Header: fully separate from illustration and labels
    rr(ax, x, y+h-.070, w, .070, fc=color, ec=color, lw=0, radius=.020, z=3)
    ax.add_patch(Rectangle((x, y+h-.050), w, .050, fc=color, ec="none", zorder=3))
    ax.add_patch(Circle((x+.027, y+h-.035), .016, fc=WHITE, ec="none", zorder=4))
    txt(ax, x+.027, y+h-.035, str(n), 8.0, color, "bold", "center", "center", 5)
    txt(ax, x+.052, y+h-.035, title, 9.2, WHITE, "bold", "left", "center", 5)
    return {
        "icon_bottom": y+h-.270,
        "icon_top": y+h-.095,
        "text_top": y+h-.295,
        "bottom": y+.020,
    }


def dna(ax, cx, cy, s=1, c1=BLUE, c2="#77C9E7"):
    ys = [cy-.060*s+i*.010*s for i in range(13)]
    a, b = [], []
    for i, yy in enumerate(ys):
        phase=i/12*2.2*math.pi
        a.append((cx-.020*s*math.sin(phase), yy))
        b.append((cx+.020*s*math.sin(phase), yy))
    ax.plot(*zip(*a), color=c1, lw=1.8, zorder=6)
    ax.plot(*zip(*b), color=c2, lw=1.8, zorder=6)
    for i in range(0,13,2):
        ax.plot([a[i][0],b[i][0]],[ys[i],ys[i]],color="#8BA9B9",lw=.7,zorder=5)


def immune_cells(ax, cx, cy):
    for dx,dy,r,fc,nc in [(-.022,.008,.025,"#77CCD5",VIOLET),(.024,-.022,.022,"#F1BC63","#8E78B8"),(.030,.030,.019,"#EEA4B4","#9774B5")]:
        ax.add_patch(Circle((cx+dx,cy+dy),r,fc=fc,ec=WHITE,lw=1.1,zorder=6))
        ax.add_patch(Circle((cx+dx-r*.12,cy+dy+r*.05),r*.42,fc=nc,ec="none",zorder=7))
        ax.add_patch(Circle((cx+dx+r*.48,cy+dy+r*.30),r*.10,fc=WHITE,ec="none",zorder=7))


def ld_and_locus(ax, cx, cy):
    n=5; step=.016
    for i in range(n):
        for j in range(i+1):
            shades=["#EDF6FA","#C9E7F2","#93D0E5","#58B2D3","#218ABA"]
            ax.add_patch(Rectangle((cx-.082+j*step,cy+.035-i*step),step*.92,step*.92,
                                   fc=shades[max(0,4-abs(i-j))],ec=WHITE,lw=.3,zorder=6))
    ax.plot([cx-.082,cx-.002],[cy+.051,cy-.029],color="#B2C3CD",lw=.8,zorder=7)
    xs=[cx+.025+i*.015 for i in range(7)]
    hs=[.012,.022,.038,.083,.048,.029,.017]
    for xx,hh in zip(xs,hs):
        ax.plot([xx,xx],[cy-.030,cy-.030+hh],color=TEAL,lw=1,zorder=6)
        ax.add_patch(Circle((xx,cy-.030+hh),.0037,fc=TEAL,ec="none",zorder=7))
    ax.add_patch(Circle((cx+.070,cy+.045),.025,fc="none",ec=AMBER,lw=1.8,zorder=8))


def databases(ax, cx, cy):
    for dx, color, label in [(-.045,AMBER,"eQTL"),(.045,TEAL,"GWAS")]:
        x=cx+dx-.025; y=cy-.035; w=.050; h=.070
        ax.add_patch(Rectangle((x,y),w,h,fc="#F2F8FA",ec=color,lw=1,zorder=5))
        ax.add_patch(Ellipse((x+w/2,y+h),w,.016,fc="#E8F2F5",ec=color,lw=1,zorder=6))
        ax.add_patch(Ellipse((x+w/2,y),w,.016,fc="#F2F8FA",ec=color,lw=1,zorder=6))
        ax.add_patch(Ellipse((x+w/2,y+h*.48),w,.014,fc="none",ec=color,lw=.7,zorder=6))
        txt(ax,cx+dx,cy+.062,label,7.1,color,"bold","center","center")


def liver_with_omics(ax, cx, cy):
    verts=[(cx-.078,cy-.025),(cx-.064,cy+.035),(cx-.010,cy+.055),(cx+.054,cy+.044),
           (cx+.073,cy+.004),(cx+.032,cy-.032),(cx-.033,cy-.040),(cx-.078,cy-.025)]
    codes=[MplPath.MOVETO]+[MplPath.CURVE3]*6+[MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts,codes),fc="#DA8262",ec="#B85F44",lw=1,zorder=5))
    ax.plot([cx+.010,cx+.008],[cy+.043,cy-.020],color="#F7C5B2",lw=1.8,zorder=6)
    for dx,dy,fc in [(-.035,-.012,"#72C7D4"),(.028,.005,"#EAA0B2")]:
        ax.add_patch(Circle((cx+dx,cy+dy),.014,fc=fc,ec=WHITE,lw=1,zorder=7))
        ax.add_patch(Circle((cx+dx-.002,cy+dy),.006,fc=VIOLET,ec="none",zorder=8))
    # RNA and protein, spatially above/right of liver
    pts=[]
    for i in range(28):
        xx=cx+.065+i*.0023; yy=cy+.028+.008*math.sin(i/3)
        pts.append((xx,yy))
    ax.plot(*zip(*pts),color=TEAL,lw=1.6,zorder=7)
    nodes=[(cx+.075,cy-.020),(cx+.095,cy-.004),(cx+.117,cy-.021),(cx+.137,cy-.002)]
    ax.plot(*zip(*nodes),color=VIOLET,lw=1.5,zorder=7)
    for p in nodes: ax.add_patch(Circle(p,.005,fc="#B7ABE0",ec=VIOLET,lw=.6,zorder=8))


def receptor_drug(ax, cx, cy):
    # receptor above TYK2, drug to right; labels stay outside icon zone
    ax.plot([cx-.026,cx+.026],[cy+.020,cy+.020],color="#7E9AAA",lw=4,zorder=5)
    for dx,c in [(-.026,BLUE),(.026,"#60BCE5")]:
        ax.plot([cx+dx,cx+dx],[cy+.020,cy+.072],color=c,lw=2.5,zorder=6)
        ax.add_patch(Circle((cx+dx,cy+.083),.012,fc="#CDEAF5",ec=c,lw=1,zorder=7))
    ax.add_patch(Ellipse((cx,cy-.024),.070,.045,fc="#EEE9F8",ec=VIOLET,lw=1.2,zorder=6))
    txt(ax,cx,cy-.024,"TYK2",5.7,VIOLET,"bold","center","center")
    ax.add_patch(FancyArrowPatch((cx+.048,cy-.002),(cx+.105,cy-.002),arrowstyle="-|>",
                                 mutation_scale=12,color="#91A2AE",lw=1.5,zorder=7))
    # capsule
    x=cx+.123; y=cy-.002
    ax.add_patch(FancyBboxPatch((x-.031,y-.010),.062,.020,boxstyle="round,pad=0,rounding_size=.010",
                                fc="#F8B49C",ec=ORANGE,lw=1.1,zorder=7))


def main():
    fig,ax=plt.subplots(figsize=(16,7.8),facecolor=BG)
    ax.set_facecolor(BG); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")

    txt(ax,.028,.958,"Single-cell eQTL-guided translational prioritization in PBC",19,INK,"bold")
    txt(ax,.029,.902,"A signal-aware evidence funnel from immune-cell regulation to a bounded therapeutic hypothesis",9.5,MUTED)
    badge(ax,.812,.907,.160,.035,"GENETICS  →  CONTEXT  →  TRANSLATION","#E6F2F8",BLUE,6.7)

    y=.264; h=.545; gap=.014
    widths=[.166,.166,.166,.190,.225]
    xs=[.027]
    for i in range(1,5): xs.append(xs[-1]+widths[i-1]+gap)
    colors=[BLUE,TEAL,AMBER,VIOLET,ORANGE]
    titles=["Discovery","Signal resolution","Replication","Tissue context","Translation gate"]
    regions=[card(ax,x,y,w,h,c,i+1,t) for i,(x,w,c,t) in enumerate(zip(xs,widths,colors,titles))]

    # 1 Discovery: strict three-zone layout
    x,w=xs[0],widths[0]; cx=x+w/2
    dna(ax,cx-.040,.658,.85); immune_cells(ax,cx+.045,.655)
    txt(ax,cx,.553,"PBC GWAS  ×  OneK1K sc-eQTL",7.8,BLUE,"bold","center")
    txt(ax,cx,.517,"56 loci  •  17 immune-cell populations",7.0,INK,"normal","center")
    badge(ax,x+.016,.302,w-.032,.038,"16 overlap loci  |  152 records","#E5F2F9",BLUE,6.4)

    # 2 Signal resolution
    x,w=xs[1],widths[1]; cx=x+w/2
    ld_and_locus(ax,cx,.646)
    txt(ax,cx,.553,"Signed LD  •  SuSiE fine-mapping",7.45,TEAL,"bold","center")
    txt(ax,cx,.516,"Signal-level coloc.susie",7.0,INK,"normal","center")
    badge(ax,x+.016,.302,w-.032,.038,"10 loci  |  23 formal tests","#E4F4F0",TEAL,6.4)

    # 3 Replication
    x,w=xs[2],widths[2]; cx=x+w/2
    databases(ax,cx,.648)
    txt(ax,cx,.553,"Molecular replication",7.35,AMBER,"bold","center")
    txt(ax,cx,.520,"Perez NK  •  GTEx blood  •  Nathan T",6.55,INK,"normal","center")
    txt(ax,cx,.477,"Disease replication",7.35,TEAL,"bold","center")
    txt(ax,cx,.444,"FinnGen  •  UK Biobank",6.75,INK,"normal","center")
    badge(ax,x+.016,.302,w-.032,.038,"Independent evidence only","#FFF2DA","#9A650E",6.2)

    # 4 Tissue context
    x,w=xs[3],widths[3]; cx=x+w/2
    liver_with_omics(ax,cx-.025,.653)
    txt(ax,cx,.553,"Liver donor-level analysis",7.45,VIOLET,"bold","center")
    txt(ax,cx,.519,"HRA008003: 5 PBC + 5 controls",6.7,INK,"normal","center")
    txt(ax,cx,.478,"External liver scRNA  •  bulk RNA",6.55,INK,"normal","center")
    txt(ax,cx,.446,"PXD083842 proteomic audit",6.55,INK,"normal","center")
    badge(ax,x+.016,.302,w-.032,.038,"Mechanistic context—not a biomarker","#EFEAF8",VIOLET,5.95)

    # 5 Translation gate
    x,w=xs[4],widths[4]; cx=x+.065
    receptor_drug(ax,cx,.647)
    txt(ax,cx,.757,"IL12RB2",7.3,BLUE,"bold","center")
    txt(ax,x+.165,.689,"Deucravacitinib",7.2,ORANGE,"bold","center")
    txt(ax,x+.018,.553,"Direction anchor",6.9,VIOLET,"bold")
    txt(ax,x+.018,.521,"TYK2 rs34536443",6.65,INK)
    badge(ax,x+.117,.517,w-.133,.045,"CONDITIONAL GO","#E3F4E9",GREEN,6.6)
    badge(ax,x+.018,.452,w-.034,.042,"Direct IL12RB2 docking: NO-GO","#FCE8E5",RED,6.25)
    txt(ax,cx+w*.22,.410,"No direct PBC efficacy trial; no validated\nIL12RB2 small-molecule pocket/activity series",6.15,MUTED,"normal","center")
    badge(ax,x+.018,.302,w-.034,.038,"Research priority—not treatment advice","#FFF0E8",ORANGE,6.05)

    # Arrows live only in gutters and cannot cover card content
    for i in range(4):
        arrow(ax,xs[i]+widths[i]+.002,xs[i+1]-.002,.579,colors[i])

    # Boundary/caveat strip: no separate decorative icons
    rr(ax,.027,.075,.945,.125,fc=WHITE,ec=LINE,lw=.9,radius=.018,z=2)
    rr(ax,.027,.075,.168,.125,fc="#344C61",ec="#344C61",lw=0,radius=.018,z=3)
    ax.add_patch(Rectangle((.175,.075),.020,.125,fc="#344C61",ec="none",zorder=3))
    txt(ax,.047,.165,"EVIDENCE BOUNDARIES",7.6,WHITE,"bold")
    txt(ax,.047,.130,"Negative and discordant results\nremain visible in the evidence chain",6.8,"#DFE9F0")
    # Three evenly spaced textual boundaries
    items=[
        ("Liver RNA","No transcriptome-wide FDR-significant\ncase–control contrast"),
        ("Protein / trial","Compartment-specific discordance;\nustekinumab 0/20 responders"),
        ("Interpretation","Counts are screened combinations,\nnot independent discoveries"),
    ]
    starts=[.222,.482,.742]
    for (head,body),sx in zip(items,starts):
        ax.add_patch(Circle((sx,.151),.006,fc={"Liver RNA":AMBER,"Protein / trial":RED,"Interpretation":BLUE}[head],ec="none",zorder=5))
        txt(ax,sx+.014,.166,head,7.0,INK,"bold")
        txt(ax,sx+.014,.135,body,6.55,MUTED)

    txt(ax,.972,.040,"CS, credible set; eQTL, expression quantitative trait locus; LD, linkage disequilibrium.",
        6.1,"#84929E","normal","right","bottom")

    fig.subplots_adjust(0,0,1,1)
    for stem in ("Fig01_study_overview_final",):
        fig.savefig(OUT/f"{stem}.png",dpi=300,facecolor=BG,bbox_inches="tight",pad_inches=.04)
        fig.savefig(OUT/f"{stem}.pdf",facecolor=BG,bbox_inches="tight",pad_inches=.04)
        fig.savefig(OUT/f"{stem}.svg",facecolor=BG,bbox_inches="tight",pad_inches=.04)
    plt.close(fig)


if __name__=="__main__":
    main()
