from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle


REPO = Path(__file__).resolve().parents[2]
FIGURE_DIR = REPO / "results" / "figures" / "supplementary"
SOURCE_DIR = REPO / "data" / "derived" / "figures" / "supplementary"

REGULATORY = SOURCE_DIR / "FigS04A_credible_set_regulatory_annotation.tsv"
GENE_LEVEL = SOURCE_DIR / "FigS04B_IL12RB2_response_state_gene_level.tsv"
ALLELE_EFFECTS = SOURCE_DIR / "FigS04C_rs11209064_activation_effects.tsv"
LD_LONG = SOURCE_DIR / "FigS04D_response_eQTL_PBC_LD.tsv"
CROSS_TRAIT = SOURCE_DIR / "FigS04E_cross_autoimmune_context.tsv"

INK = "#25364A"
MUTED = "#6E7B87"
GRID = "#DCE4EA"
LIGHT = "#F5F7F9"
BLUE = "#2878B5"
TEAL = "#2A9D8F"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
PURPLE = "#7A5AA6"
L1 = "#4C78A8"
L2 = "#E45756"
L3 = "#72B7B2"


def panel_label(ax, letter):
    ax.text(-0.11, 1.08, letter, transform=ax.transAxes, fontsize=13.5,
            fontweight="bold", color=INK, va="top")


def clean_axes(ax, grid_axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=INK, labelsize=8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.7, alpha=0.75, zorder=0)


def panel_a(ax, regulatory):
    data = regulatory.copy()
    order = ["rs6679356", "rs6702599", "rs3828069", "rs3790569", "rs11209064", "rs112588573", "rs72932868"]
    data["rsid"] = pd.Categorical(data["rsid"], categories=order, ordered=True)
    data = data.sort_values("rsid").reset_index(drop=True)
    component_color = {"L1": L1, "L2": L2, "L3": L3}
    y_positions = np.arange(len(data))[::-1]
    variant_x = 0.05
    component_x = 1.55
    pip_left = 2.15
    pip_width = 0.72
    pip_value_x = 3.05
    ccre_x = 4.02
    e2g_x = 5.18

    ax.axhspan(6.48, 7.28, color="#EEF2F5", zorder=0)
    for xpos in [1.25, 1.92, 3.42, 4.60]:
        ax.axvline(xpos, color="white", linewidth=1.4, zorder=1)

    headers = [
        (variant_x, "Variant", "left"),
        (component_x, "Signal", "center"),
        ((pip_left + pip_value_x) / 2, "PBC PIP", "center"),
        (ccre_x, "ENCODE\ncCRE", "center"),
        (e2g_x, "Activated T-cell\nIL12RB2 E2G", "center"),
    ]
    for xpos, label, alignment in headers:
        ax.text(xpos, 6.88, label, ha=alignment, va="center", fontsize=8.4,
                fontweight="bold", color=INK, linespacing=1.0)

    for y, row in zip(y_positions, data.itertuples()):
        color = component_color[row.component]
        pip = float(row.pbc_pip)
        ax.axhspan(y - 0.44, y + 0.44, color=color, alpha=0.065, zorder=0)
        ax.text(variant_x, y, str(row.rsid), ha="left", va="center", fontsize=8.8,
                color=INK, fontweight="bold")

        signal_box = FancyBboxPatch(
            (component_x - 0.22, y - 0.23), 0.44, 0.46,
            boxstyle="round,pad=0.02,rounding_size=0.07",
            linewidth=0, facecolor=color, zorder=3,
        )
        ax.add_patch(signal_box)
        ax.text(component_x, y, str(row.component), ha="center", va="center",
                fontsize=8.3, color="white", fontweight="bold", zorder=4)

        ax.add_patch(Rectangle((pip_left, y - 0.13), pip_width, 0.26,
                               linewidth=0, facecolor="#DCE4EA", zorder=2))
        ax.add_patch(Rectangle((pip_left, y - 0.13), max(0.025, pip_width * pip), 0.26,
                               linewidth=0, facecolor=color, zorder=3))
        ax.text(pip_value_x, y, f"{pip:.2f}", ha="center", va="center",
                fontsize=8.5, color=INK, fontweight="bold")

        if int(row.encode_ccre_count) > 0:
            ccre_box = FancyBboxPatch(
                (ccre_x - 0.33, y - 0.22), 0.66, 0.44,
                boxstyle="round,pad=0.02,rounding_size=0.06",
                linewidth=0, facecolor=TEAL, zorder=3,
            )
            ax.add_patch(ccre_box)
            ax.text(ccre_x, y, "dELS", ha="center", va="center", fontsize=8.0,
                    color="white", fontweight="bold", zorder=4)
        else:
            ax.text(ccre_x, y, "—", ha="center", va="center", fontsize=9.2, color=MUTED)

        count = int(row.e2g_il12rb2_t_cell_count)
        if count:
            e2g_box = FancyBboxPatch(
                (e2g_x - 0.25, y - 0.24), 0.50, 0.48,
                boxstyle="round,pad=0.02,rounding_size=0.12",
                linewidth=0, facecolor=ORANGE, zorder=3,
            )
            ax.add_patch(e2g_box)
            ax.text(e2g_x, y, str(count), ha="center", va="center", fontsize=8.8,
                    color="white", fontweight="bold", zorder=4)
        else:
            ax.text(e2g_x, y, "—", ha="center", va="center", fontsize=9.2, color=MUTED)

    ax.set_xlim(-0.08, 5.72)
    ax.set_ylim(-0.52, 7.30)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Credible-set regulatory architecture", loc="left", fontweight="bold",
                 fontsize=10.7, pad=12)
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_label(ax, "A")


def classify_state(row):
    if row.study_label == "Randolph_2021":
        return "Influenza / naive"
    condition = str(row.condition_label).lower()
    if condition == "naive":
        return "Unstimulated CD4"
    if "16h" in condition:
        return "anti-CD3/CD28 16 h"
    if "40h" in condition:
        return "anti-CD3/CD28 40 h"
    if "5d" in condition:
        return "anti-CD3/CD28 5 d"
    return "Other"


def short_state_label(sample_group):
    replacements = {
        "combined_CD4_Memory_STIM_40H": "Memory CD4, 40 h",
        "combined_CD4_Memory_STIM_5D": "Memory CD4, 5 d",
        "TEM_STIM_5D": "TEM, 5 d",
        "TCM_STIM_40H": "TCM, 40 h",
        "TN_A_STIM_16H": "TN-A, 16 h",
    }
    return replacements.get(sample_group, sample_group)


def panel_b(ax, gene_level):
    data = gene_level.copy()
    data["state_class"] = data.apply(classify_state, axis=1)
    data["minus_log10_p_beta"] = -np.log10(data["p_beta"].astype(float).clip(lower=1e-300))
    categories = ["Influenza / naive", "Unstimulated CD4", "anti-CD3/CD28 16 h", "anti-CD3/CD28 40 h", "anti-CD3/CD28 5 d"]
    palette = {
        "Influenza / naive": "#A9B3BC",
        "Unstimulated CD4": "#7E8A95",
        "anti-CD3/CD28 16 h": "#56B4E9",
        "anti-CD3/CD28 40 h": VERMILION,
        "anti-CD3/CD28 5 d": ORANGE,
    }
    rng = np.random.default_rng(20261026)
    for xpos, category in enumerate(categories):
        subset = data[data.state_class == category].copy()
        jitter = rng.uniform(-0.20, 0.20, len(subset))
        if category == "anti-CD3/CD28 5 d":
            fixed_offsets = {
                "combined_CD4_Memory_STIM_5D": -0.10,
                "TEM_STIM_5D": 0.10,
            }
            for sample_group, offset in fixed_offsets.items():
                jitter[subset.sample_group.eq(sample_group).to_numpy()] = offset
        x_values = np.full(len(subset), xpos) + jitter
        positive = subset.p_beta.astype(float) < 0.05
        ax.scatter(x_values, subset.minus_log10_p_beta,
                   s=np.where(positive, 55, 25), color=np.where(positive, VERMILION, palette[category]),
                   alpha=np.where(positive, 1.0, 0.72), edgecolor="white", linewidth=0.6, zorder=3)
        label_offsets = {
            "combined_CD4_Memory_STIM_40H": (6, 7, "left"),
            "TCM_STIM_40H": (7, -14, "left"),
            "TN_A_STIM_16H": (6, 6, "left"),
            "combined_CD4_Memory_STIM_5D": (-10, 13, "right"),
            "TEM_STIM_5D": (10, 4, "left"),
        }
        for x, (_, row) in zip(x_values, subset.iterrows()):
            if float(row.p_beta) < 0.05:
                dx, dy, alignment = label_offsets[row.sample_group]
                leader = None
                if row.sample_group in {"combined_CD4_Memory_STIM_5D", "TEM_STIM_5D"}:
                    leader = dict(arrowstyle="-", color=MUTED, linewidth=0.65,
                                  shrinkA=1.5, shrinkB=4.0)
                ax.annotate(short_state_label(row.sample_group), (x, row.minus_log10_p_beta),
                            xytext=(dx, dy), textcoords="offset points", fontsize=6.6,
                            color=INK, rotation=22, ha=alignment, va="bottom",
                            arrowprops=leader, zorder=5)
    threshold = -np.log10(0.05)
    ax.axhline(threshold, color=INK, linestyle="--", linewidth=1.0)
    ax.text(0.02, threshold + 0.10, r"$p_{\beta}=0.05$", transform=ax.get_yaxis_transform(),
            fontsize=7.2, color=INK, ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.94))
    ax.set_xticks(range(len(categories)), ["Influenza /\nnaive", "Unstimulated\nCD4", "16 h", "40 h", "5 d"])
    ax.set_ylabel("−log10(gene-level p_beta)")
    ax.set_title("IL12RB2 across 56 immune states", loc="left", fontweight="bold", fontsize=10.5)
    ax.text(0.02, 0.97, "5/56 states with p_beta<0.05\nall anti-CD3/CD28-activated CD4 states",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.3, color=INK,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=GRID, alpha=0.95))
    clean_axes(ax)
    panel_label(ax, "B")


def panel_c(ax, allele_effects):
    keep = allele_effects[
        (allele_effects.pbc_rsid == "rs11209064")
        & allele_effects.sample_group.isin([
            "combined_CD4_Memory_UNS_16H",
            "combined_CD4_Memory_STIM_40H",
            "TCM_UNS_16H",
            "TCM_STIM_40H",
        ])
    ].copy()
    label_map = {
        "combined_CD4_Memory_UNS_16H": "Memory CD4, unstimulated",
        "combined_CD4_Memory_STIM_40H": "Memory CD4, 40 h activated",
        "TCM_UNS_16H": "TCM, unstimulated",
        "TCM_STIM_40H": "TCM, 40 h activated",
    }
    order = list(label_map)
    keep["sample_group"] = pd.Categorical(keep.sample_group, categories=order, ordered=True)
    keep = keep.sort_values("sample_group")
    keep["label"] = keep.sample_group.map(label_map)
    keep["beta_risk"] = keep.IL12RB2_beta_per_pbc_risk_increasing_allele.astype(float)
    keep["se"] = keep.eqtl_se.astype(float)
    keep["lower"] = keep.beta_risk - 1.96 * keep.se
    keep["upper"] = keep.beta_risk + 1.96 * keep.se
    y = np.arange(len(keep))[::-1]
    activated = keep.condition_label.str.contains("40h", case=False)
    colors = np.where(activated, VERMILION, "#8996A3")
    for ypos, (_, row), color in zip(y, keep.iterrows(), colors):
        ax.errorbar(row.beta_risk, ypos,
                    xerr=[[row.beta_risk - row.lower], [row.upper - row.beta_risk]],
                    fmt="none", ecolor=color, elinewidth=1.8, capsize=3, zorder=2)
        ax.scatter(row.beta_risk, ypos, s=55, color=color, edgecolor="white", linewidth=0.7, zorder=3)
    ax.axvline(0, color=INK, linewidth=0.9, linestyle="--")
    for ypos, row in zip(y, keep.itertuples()):
        ax.text(1.25, ypos, f"P={float(row.eqtl_pvalue):.2g}", fontsize=6.8, color=INK, va="center", ha="right")
    ax.set_yticks(y, keep.label)
    ax.set_xlabel("IL12RB2 eQTL effect per PBC risk allele G")
    ax.set_xlim(-0.65, 1.34)
    ax.set_title("rs11209064 direction by activation state", loc="left", fontweight="bold", fontsize=10.5)
    clean_axes(ax, grid_axis="x")
    panel_label(ax, "C")


def panel_d(ax, ld_long):
    row_order = ["rs2201584", "rs3790561"]
    col_order = ["rs6679356", "rs6702599", "rs3828069", "rs3790569", "rs11209064", "rs112588573", "rs72932868"]
    matrix = ld_long.pivot(index="response_susie_variant", columns="pbc_credible_variant", values="signed_r_ref_based")
    matrix = matrix.reindex(index=row_order, columns=col_order).astype(float)
    values = matrix.to_numpy()
    im = ax.imshow(values, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), aspect="auto")
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            value = values[row, col]
            ax.text(col, row, f"{value:.2f}\nr²={value**2:.2f}", ha="center", va="center",
                    fontsize=6.6, color="white" if abs(value) > 0.62 else INK,
                    fontweight="bold" if abs(value) > 0.8 else "normal")
    ax.set_xticks(np.arange(len(col_order)), col_order, rotation=38, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(len(row_order)), row_order)
    ax.tick_params(length=0)
    ax.set_ylim(1.5, -0.90)
    group_headers = [
        (-0.5, 2.0, 0.5, "L1", L1),
        (1.5, 3.0, 3.0, "L2", L2),
        (4.5, 2.0, 5.5, "L3", L3),
    ]
    for left, width, center, label, color in group_headers:
        ax.add_patch(Rectangle((left, -0.88), width, 0.34,
                               facecolor=color, alpha=0.12, edgecolor="white",
                               linewidth=0.8, clip_on=True, zorder=4))
        ax.text(center, -0.71, label, color=color, fontsize=8.0,
                fontweight="bold", ha="center", va="center", zorder=5)
    ax.set_title("Response-eQTL credible variants\nversus PBC components", loc="left",
                 fontweight="bold", fontsize=10.5, pad=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Signed r", fontsize=7.5, color=INK)
    cbar.ax.tick_params(labelsize=7, colors=INK)
    panel_label(ax, "D")


def format_p(value):
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def panel_e(ax, cross_trait):
    data = cross_trait.set_index("trait").loc[["PSC", "RA", "IBD", "Psoriasis"]].reset_index()
    data["max_r2"] = data[["lead_r2_to_PBC_L1", "lead_r2_to_PBC_L2", "lead_r2_to_PBC_L3"]].max(axis=1, skipna=True)
    column_lefts = np.array([0.00, 1.10, 2.70, 4.18, 5.42, 6.66])
    column_widths = np.array([1.10, 1.60, 1.48, 1.24, 1.24, 3.24])
    headers = ["Trait", "Regional lead P", "max r² to PBC", "PP.H3", "PP.H4", "Signal interpretation"]
    row_y = np.array([3.45, 2.45, 1.45, 0.45])
    row_height = 0.82
    header_y = 4.35

    for left, width, header in zip(column_lefts, column_widths, headers):
        ax.add_patch(Rectangle((left, header_y - 0.36), width, 0.72,
                               facecolor="#EAF0F4", edgecolor="white", linewidth=1.1))
        ax.text(left + 0.10 if header == "Trait" else left + width / 2, header_y, header,
                ha="left" if header == "Trait" else "center", va="center",
                fontsize=8.9, fontweight="bold", color=INK)

    for index, (ypos, row) in enumerate(zip(row_y, data.itertuples())):
        base_color = "#F8FAFB" if index % 2 == 0 else "#F1F5F7"
        for left, width in zip(column_lefts, column_widths):
            ax.add_patch(Rectangle((left, ypos - row_height / 2), width, row_height,
                                   facecolor=base_color, edgecolor="white", linewidth=1.1, zorder=0))

        ax.text(column_lefts[0] + 0.10, ypos, row.trait, ha="left", va="center",
                fontsize=9.3, fontweight="bold", color=INK)

        pvalue = float(row.regional_lead_pvalue)
        p_text = format_p(pvalue)
        p_fill = "#FCE8DD" if bool(row.genome_wide_significant_region) else "#E9EEF2"
        p_color = VERMILION if bool(row.genome_wide_significant_region) else INK
        ax.add_patch(Rectangle((column_lefts[1], ypos - row_height / 2), column_widths[1], row_height,
                               facecolor=p_fill, edgecolor="white", linewidth=1.1, zorder=1))
        ax.text(column_lefts[1] + column_widths[1] / 2, ypos, f"P = {p_text}",
                ha="center", va="center", fontsize=8.8, color=p_color,
                fontweight="bold" if bool(row.genome_wide_significant_region) else "normal", zorder=2)

        max_r2 = float(row.max_r2) if pd.notna(row.max_r2) else np.nan
        r2_text = f"{max_r2:.3f}" if np.isfinite(max_r2) else "Not available"
        ax.add_patch(Rectangle((column_lefts[2], ypos - row_height / 2), column_widths[2], row_height,
                               facecolor="#E4F0F7" if np.isfinite(max_r2) else "#F1F3F5",
                               edgecolor="white", linewidth=1.1, zorder=1))
        ax.text(column_lefts[2] + column_widths[2] / 2, ypos, r2_text,
                ha="center", va="center", fontsize=8.8,
                color=BLUE if np.isfinite(max_r2) else MUTED, fontweight="bold", zorder=2)

        for column_index, value, color in [
            (3, float(row.coloc_PP_H3_default), PURPLE),
            (4, float(row.coloc_PP_H4_default), TEAL),
        ]:
            ax.add_patch(Rectangle(
                (column_lefts[column_index], ypos - row_height / 2),
                column_widths[column_index], row_height,
                facecolor=color, alpha=0.16 + 0.80 * value,
                edgecolor="white", linewidth=1.1, zorder=1,
            ))
            ax.text(column_lefts[column_index] + column_widths[column_index] / 2, ypos,
                    f"{value:.3f}", ha="center", va="center", fontsize=9.0,
                    color="white" if value > 0.55 else INK,
                    fontweight="bold", zorder=2)

        is_distinct = row.classification == "distinct_local_component"
        interpretation = "Distinct IL23R-side\ncomponent" if is_distinct else "No detectable\nregional association"
        interpretation_fill = "#FCE8DD" if is_distinct else "#EEF2F5"
        interpretation_color = VERMILION if is_distinct else MUTED
        ax.add_patch(Rectangle((column_lefts[5], ypos - row_height / 2), column_widths[5], row_height,
                               facecolor=interpretation_fill, edgecolor="white", linewidth=1.1, zorder=1))
        ax.text(column_lefts[5] + 0.16, ypos, interpretation, ha="left", va="center",
                fontsize=8.8, linespacing=1.05, color=interpretation_color,
                fontweight="bold" if is_distinct else "normal", zorder=2)

    ax.set_xlim(0, 9.90)
    ax.set_ylim(-0.75, 4.77)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Cross-autoimmune context at 1p31.3", loc="left", fontweight="bold",
                 fontsize=10.7, pad=12)
    ax.text(0.0, -0.48, r"Default coloc prior $p_{12}=1\times10^{-5}$. IBD and psoriasis support H3, not a shared PBC signal.",
            fontsize=8.2, color=INK, va="center")
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_label(ax, "E")


def export_source_data(regulatory, gene_level, allele_effects, ld_long, cross_trait):
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    regulatory.to_csv(SOURCE_DIR / "FigS04A_credible_set_regulatory_annotation.tsv", sep="\t", index=False)
    gene_level.to_csv(SOURCE_DIR / "FigS04B_IL12RB2_response_state_gene_level.tsv", sep="\t", index=False)
    allele_effects[
        (allele_effects.pbc_rsid == "rs11209064")
        & allele_effects.sample_group.isin([
            "combined_CD4_Memory_UNS_16H",
            "combined_CD4_Memory_STIM_40H",
            "TCM_UNS_16H",
            "TCM_STIM_40H",
        ])
    ].to_csv(SOURCE_DIR / "FigS04C_rs11209064_activation_effects.tsv", sep="\t", index=False)
    ld_long.to_csv(SOURCE_DIR / "FigS04D_response_eQTL_PBC_LD.tsv", sep="\t", index=False)
    cross_trait.to_csv(SOURCE_DIR / "FigS04E_cross_autoimmune_context.tsv", sep="\t", index=False)


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    regulatory = pd.read_csv(REGULATORY, sep="\t")
    gene_level = pd.read_csv(GENE_LEVEL, sep="\t")
    allele_effects = pd.read_csv(ALLELE_EFFECTS, sep="\t")
    ld_long = pd.read_csv(LD_LONG, sep="\t")
    cross_trait = pd.read_csv(CROSS_TRAIT, sep="\t")
    export_source_data(regulatory, gene_level, allele_effects, ld_long, cross_trait)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8.3,
        "axes.titlesize": 10.5,
        "axes.labelsize": 8.3,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig = plt.figure(figsize=(13.2, 12.8), constrained_layout=False)
    grid = fig.add_gridspec(3, 12, left=0.060, right=0.985, top=0.905, bottom=0.06,
                            hspace=0.68, wspace=0.90, height_ratios=[1.08, 0.92, 0.88])
    ax_a = fig.add_subplot(grid[0, 0:6])
    ax_b = fig.add_subplot(grid[0, 6:12])
    ax_c = fig.add_subplot(grid[1, 0:5])
    ax_d = fig.add_subplot(grid[1, 5:12])
    ax_e = fig.add_subplot(grid[2, 0:12])

    panel_a(ax_a, regulatory)
    panel_b(ax_b, gene_level)
    panel_c(ax_c, allele_effects)
    panel_d(ax_d, ld_long)
    panel_e(ax_e, cross_trait)

    fig.suptitle("Supplementary Figure S4 | A state-resolved secondary IL12RB2 regulatory component",
                 x=0.02, y=0.985, ha="left", fontsize=14.2, fontweight="bold", color=INK)
    fig.text(0.985, 0.015,
             "Activation-state eQTL and LD compatibility are supportive; no formal genotype-by-stimulation interaction was established.",
             ha="right", fontsize=7.4, color=VERMILION)

    stem = FIGURE_DIR / "FigS04_IL12RB2_activation_and_cross_autoimmune_context_final"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
