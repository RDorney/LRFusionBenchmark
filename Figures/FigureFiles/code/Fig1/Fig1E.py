"""Figure 1E. Fusion-type recall against false-call profile.

Left panel: recovered spiked-in fusions by fusion type (intra-chromosomal,
inter-chromosomal, read-through, tri-fusion), out of 100 per type at 100Gb /
95% identity. Right panel: false calls split into chimera and other errors,
where other errors combines the mitochondrial-genomic, mitochondrial,
self-misalignment and sense-antisense classes. Independent y-axes sit on the
outer edges and the six callers share one legend.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig1'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

# ----- Data -----
df = pd.read_csv(DATA / 'fusion_profile_stat_summary.tsv', sep='\t')
df = df[(df.control == 'positive') & (df.depth == '100GB') & (df.Sequence_Identity == '95%')]

TRUE_TYPES = [
    ('intra_chromosome', 'Intra-\nchromosomal'),
    ('hybrid',           'Inter-\nchromosomal'),
    ('read_through',     'Read-through'),
    ('tri_fusion',       'Tri-fusion\n(3-gene)'),
]

# Right panel: plain chimera, and all other error subtypes merged into one.
RAW_FALSE_KEYS = {
    'chimera':      ['false_fusion'],
    'other_errors': ['false_fusion:mitochondrial_genomic',
                     'false_fusion:mitochondrial',
                     'false_fusion:self_misalignment',
                     'false_fusion:Sense-Antisense'],
}
FALSE_TYPES_MERGED = [
    ('chimera',      'Chimera'),
    ('other_errors', 'Other errors'),
]

mat_T = (df[df.recall_category == 'True_Recall']
            .pivot(index='Algorithm', columns='fusionType', values='fusion_type_count')
            .reindex(index=TOOL_ORDER, columns=[k for k, _ in TRUE_TYPES])
            .fillna(0).astype(int))

# Build collapsed false matrix
false_raw = (df[df.recall_category == 'False_Call']
                .pivot(index='Algorithm', columns='fusionType', values='fusion_type_count')
                .reindex(index=TOOL_ORDER)
                .fillna(0).astype(int))

mat_F = pd.DataFrame(index=TOOL_ORDER, dtype=int)
for bucket_key, raw_keys in RAW_FALSE_KEYS.items():
    cols = [c for c in raw_keys if c in false_raw.columns]
    mat_F[bucket_key] = false_raw[cols].sum(axis=1) if cols else 0
mat_F = mat_F[[k for k, _ in FALSE_TYPES_MERGED]]

Y_T_MAX = 55
Y_F_MAX = int(np.ceil(max(mat_F.values.max(), 1) / 100.0) * 100)
n_tools = len(TOOL_ORDER)


# ----- Figure -----
fig = plt.figure(figsize=(12.0, 5.2))
fig.patch.set_facecolor('white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_axis_off()
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

L_x0, L_x1 = 0.080, 0.480
R_x0, R_x1 = 0.500, 0.910
plot_bottom = 0.155
plot_top    = 0.760
bracket_y      = 0.798
bracket_text_y = 0.810


def yfrac_T(v):
    return plot_bottom + v / Y_T_MAX * (plot_top - plot_bottom)

def yfrac_F(v):
    return plot_bottom + v / Y_F_MAX * (plot_top - plot_bottom)


def draw_bracket(x0, x1, label):
    ax.plot([x0, x1], [bracket_y, bracket_y],
            transform=fig.transFigure,
            color=NEUTRAL['rule'], linewidth=0.9, zorder=4)
    for x in (x0, x1):
        ax.plot([x, x], [bracket_y, bracket_y - 0.012],
                transform=fig.transFigure,
                color=NEUTRAL['rule'], linewidth=0.9, zorder=4)
    ax.text((x0 + x1) / 2, bracket_text_y, label,
            transform=fig.transFigure,
            ha='center', va='bottom',
            fontsize=11.5, fontweight='semibold', color=NEUTRAL['ink'])


def draw_y_axis(side, panel_x0, panel_x1, y_max, ticks):
    for tv in ticks:
        y = plot_bottom + tv / y_max * (plot_top - plot_bottom)
        ax.plot([panel_x0, panel_x1], [y, y],
                transform=fig.transFigure,
                color=NEUTRAL['grid'], linewidth=0.5, zorder=1)
        if side == 'left':
            ax.text(panel_x0 - 0.008, y, str(tv),
                    transform=fig.transFigure,
                    ha='right', va='center',
                    fontsize=8.6, color=NEUTRAL['body'])
        else:
            ax.text(panel_x1 + 0.008, y, str(tv),
                    transform=fig.transFigure,
                    ha='left', va='center',
                    fontsize=8.6, color=NEUTRAL['body'])


def draw_grouped_bars(panel_x0, panel_x1, types_def, mat, yfrac):
    n_cat = len(types_def)
    slot_w = (panel_x1 - panel_x0) / n_cat
    inner_pad = slot_w * 0.10
    bars_total_w = slot_w - 2 * inner_pad
    bar_w = bars_total_w / n_tools

    for c_idx, (key, label) in enumerate(types_def):
        slot_l = panel_x0 + c_idx * slot_w
        group_mid = slot_l + slot_w / 2
        if c_idx % 2 == 0:
            ax.add_patch(Rectangle(
                (slot_l, plot_bottom), slot_w, plot_top - plot_bottom,
                transform=fig.transFigure,
                facecolor='#FAFAF6', edgecolor='none', zorder=0))
        ax.text(group_mid, plot_bottom - 0.018, label,
                transform=fig.transFigure,
                ha='center', va='top',
                fontsize=9.2, fontweight='semibold', color=NEUTRAL['ink'])

        for t_idx, tool in enumerate(TOOL_ORDER):
            x0 = slot_l + inner_pad + t_idx * bar_w
            x_mid = x0 + bar_w / 2
            count = int(mat.loc[tool, key])
            if count > 0:
                y1 = yfrac(count)
                ax.add_patch(Rectangle(
                    (x0, plot_bottom), bar_w * 0.92, y1 - plot_bottom,
                    transform=fig.transFigure,
                    facecolor=TOOL_COLORS[tool], edgecolor='white', linewidth=0.3,
                    zorder=3))
                ax.text(x_mid, y1 + 0.008, str(count),
                        transform=fig.transFigure,
                        ha='center', va='bottom',
                        fontsize=7.4, color=NEUTRAL['body'])
            else:
                ax.plot([x0 + 0.0008, x0 + bar_w * 0.92 - 0.0008],
                        [plot_bottom + 0.001, plot_bottom + 0.001],
                        transform=fig.transFigure,
                        color=NEUTRAL['na_text'], linewidth=0.9, zorder=3,
                        solid_capstyle='butt')


# ----- LEFT panel -----
draw_bracket(L_x0, L_x1, 'True Recall')
draw_y_axis('left',  L_x0, L_x1, Y_T_MAX, ticks=[0, 10, 20, 30, 40, 50])
draw_grouped_bars(L_x0, L_x1, TRUE_TYPES, mat_T, yfrac_T)
ax.plot([L_x0, L_x1], [plot_bottom, plot_bottom],
        transform=fig.transFigure,
        color=NEUTRAL['rule'], linewidth=0.6, zorder=4)
ax.text(L_x0 - 0.045, (plot_bottom + plot_top) / 2,
        'Count',
        transform=fig.transFigure,
        rotation=90, ha='center', va='center',
        fontsize=9.5, color=NEUTRAL['body'])


# ----- RIGHT panel -----
draw_bracket(R_x0, R_x1, 'False Calls')
draw_y_axis('right', R_x0, R_x1, Y_F_MAX, ticks=list(range(0, Y_F_MAX + 1, 100)))
draw_grouped_bars(R_x0, R_x1, FALSE_TYPES_MERGED, mat_F, yfrac_F)
ax.plot([R_x0, R_x1], [plot_bottom, plot_bottom],
        transform=fig.transFigure,
        color=NEUTRAL['rule'], linewidth=0.6, zorder=4)
ax.text(R_x1 + 0.052, (plot_bottom + plot_top) / 2,
        'Count',
        transform=fig.transFigure,
        rotation=270, ha='center', va='center',
        fontsize=9.5, color=NEUTRAL['body'])


# ----- Shared tool legend (top) -----
leg_y = 0.880
total_w = R_x1 - L_x0
swatch_w, swatch_h = 0.014, 0.018
spacing = total_w / n_tools
for i, tool in enumerate(TOOL_ORDER):
    lx = L_x0 + i * spacing + 0.012
    ax.add_patch(Rectangle((lx, leg_y - swatch_h/2), swatch_w, swatch_h,
                           transform=fig.transFigure,
                           facecolor=TOOL_COLORS[tool],
                           edgecolor='white', linewidth=0.4, zorder=3))
    ax.text(lx + swatch_w + 0.005, leg_y, tool,
            transform=fig.transFigure,
            ha='left', va='center',
            fontsize=8.7, color=NEUTRAL['body'])


# ----- Title -----
ax.text(0.040, 0.972,
        'Fusion-type recall vs false-call profile',
        transform=fig.transFigure,
        ha='left', va='top',
        fontsize=13, fontweight='semibold', color=NEUTRAL['ink'])


save_figure(fig, 'Fig1E', OUTDIR)
plt.close(fig)
