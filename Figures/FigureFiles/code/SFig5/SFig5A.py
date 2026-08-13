"""Supplementary Figure 5A. Fusion-type recall against false-call profile by error class.

A two-panel mirror layout sharing one tool legend. The left panel shows true
recall by fusion type (intra-chromosomal, inter-chromosomal, read-through,
tri-fusion) on a 0 to 55 count axis; the right panel shows false calls by error
subtype (chimera, mitochondrial-genomic, mitochondrial-mitochondrial,
self-misalignment, sense-antisense) on a 0 to 700 count axis. The y-axes sit at
the outer edges so the central gap divides the two halves, and the scales are
independent because true-recall counts cap near 50 while false calls reach
several hundred. Conditions: 100 Gb depth, 95 per cent mean identity, 100
spiked fusions of each type.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'SFig5'
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
FALSE_TYPES = [
    ('false_fusion',                       'Chimera'),
    ('false_fusion:mitochondrial_genomic', 'Mito-Genomic'),
    ('false_fusion:mitochondrial',         'Mito-Mito'),
    ('false_fusion:self_misalignment',     'Self-misalign.'),
    ('false_fusion:Sense-Antisense',       'Sense-Antisense'),
]

mat_T = (df[df.recall_category == 'True_Recall']
            .pivot(index='Algorithm', columns='fusionType', values='fusion_type_count')
            .reindex(index=TOOL_ORDER, columns=[k for k, _ in TRUE_TYPES])
            .fillna(0).astype(int))
mat_F = (df[df.recall_category == 'False_Call']
            .pivot(index='Algorithm', columns='fusionType', values='fusion_type_count')
            .reindex(index=TOOL_ORDER, columns=[k for k, _ in FALSE_TYPES])
            .fillna(0).astype(int))

Y_T_MAX = 55
Y_F_MAX = 700
n_tools = len(TOOL_ORDER)


# ----- Figure -----
fig = plt.figure(figsize=(12.0, 5.2))
fig.patch.set_facecolor('white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_axis_off()
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

# Panel boundaries (figure fraction) — panels nearly touch; the bracket
# headers above each panel and the outer y-axes mark the divide.
L_x0, L_x1 = 0.080, 0.480
R_x0, R_x1 = 0.500, 0.910
plot_bottom = 0.155
plot_top    = 0.760

# header bracket y positions
bracket_y      = 0.798
bracket_text_y = 0.810


def yfrac_T(v):
    return plot_bottom + v / Y_T_MAX * (plot_top - plot_bottom)

def yfrac_F(v):
    return plot_bottom + v / Y_F_MAX * (plot_top - plot_bottom)


# ----- helper: draw bracket header above a panel -----
def draw_bracket(x0, x1, label):
    # bracket: horizontal line with downward end ticks
    ax.plot([x0, x1], [bracket_y, bracket_y],
            transform=fig.transFigure,
            color=NEUTRAL['rule'], linewidth=0.9, zorder=4)
    tick_drop = 0.012
    for x in (x0, x1):
        ax.plot([x, x], [bracket_y, bracket_y - tick_drop],
                transform=fig.transFigure,
                color=NEUTRAL['rule'], linewidth=0.9, zorder=4)
    ax.text((x0 + x1) / 2, bracket_text_y, label,
            transform=fig.transFigure,
            ha='center', va='bottom',
            fontsize=11.5, fontweight='semibold', color=NEUTRAL['ink'])


# ----- helper: draw y-axis ticks/grid for a panel -----
def draw_y_axis(side, panel_x0, panel_x1, y_max, ticks):
    """side: 'left' = ticks/labels on left of panel_x0; 'right' = on right of panel_x1."""
    for tv in ticks:
        if side == 'left':
            y = plot_bottom + tv / y_max * (plot_top - plot_bottom)
            # gridline across panel
            ax.plot([panel_x0, panel_x1], [y, y],
                    transform=fig.transFigure,
                    color=NEUTRAL['grid'], linewidth=0.5, zorder=1)
            ax.text(panel_x0 - 0.008, y, str(tv),
                    transform=fig.transFigure,
                    ha='right', va='center',
                    fontsize=8.6, color=NEUTRAL['body'])
        else:
            y = plot_bottom + tv / y_max * (plot_top - plot_bottom)
            ax.plot([panel_x0, panel_x1], [y, y],
                    transform=fig.transFigure,
                    color=NEUTRAL['grid'], linewidth=0.5, zorder=1)
            ax.text(panel_x1 + 0.008, y, str(tv),
                    transform=fig.transFigure,
                    ha='left', va='center',
                    fontsize=8.6, color=NEUTRAL['body'])


# ----- helper: draw grouped bars in a panel -----
def draw_grouped_bars(panel_x0, panel_x1, types_def, mat, yfrac):
    """Draw 6 dodged bars per category across the panel."""
    n_cat = len(types_def)
    slot_w = (panel_x1 - panel_x0) / n_cat
    inner_pad = slot_w * 0.10
    bars_total_w = slot_w - 2 * inner_pad
    bar_w = bars_total_w / n_tools

    for c_idx, (key, label) in enumerate(types_def):
        slot_l = panel_x0 + c_idx * slot_w
        group_mid = slot_l + slot_w / 2
        # alternate group banding (subtle)
        if c_idx % 2 == 0:
            ax.add_patch(Rectangle(
                (slot_l, plot_bottom), slot_w, plot_top - plot_bottom,
                transform=fig.transFigure,
                facecolor='#FAFAF6', edgecolor='none', zorder=0))
        # x-axis category label (bold so it reads as a section name, not a number)
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
                    (x0, plot_bottom), bar_w * 0.90, y1 - plot_bottom,
                    transform=fig.transFigure,
                    facecolor=TOOL_COLORS[tool], edgecolor='white', linewidth=0.25,
                    zorder=3))
                ax.text(x_mid, y1 + 0.008, str(count),
                        transform=fig.transFigure,
                        ha='center', va='bottom',
                        fontsize=7.0, color=NEUTRAL['body'])
            else:
                # zero indicator
                ax.plot([x0 + 0.0008, x0 + bar_w * 0.90 - 0.0008],
                        [plot_bottom + 0.0010, plot_bottom + 0.0010],
                        transform=fig.transFigure,
                        color=NEUTRAL['na_text'], linewidth=0.9, zorder=3,
                        solid_capstyle='butt')


# ----- LEFT panel -----
draw_bracket(L_x0, L_x1, 'True Recall')
draw_y_axis('left',  L_x0, L_x1, Y_T_MAX, ticks=[0, 10, 20, 30, 40, 50])
draw_grouped_bars(L_x0, L_x1, TRUE_TYPES, mat_T, yfrac_T)
# baseline rule
ax.plot([L_x0, L_x1], [plot_bottom, plot_bottom],
        transform=fig.transFigure,
        color=NEUTRAL['rule'], linewidth=0.6, zorder=4)
# y-axis caption (rotated, far left)
ax.text(L_x0 - 0.045, (plot_bottom + plot_top) / 2,
        'Count',
        transform=fig.transFigure,
        rotation=90, ha='center', va='center',
        fontsize=9.5, color=NEUTRAL['body'])


# ----- RIGHT panel -----
draw_bracket(R_x0, R_x1, 'False Calls')
draw_y_axis('right', R_x0, R_x1, Y_F_MAX, ticks=[0, 100, 200, 300, 400, 500, 600, 700])
draw_grouped_bars(R_x0, R_x1, FALSE_TYPES, mat_F, yfrac_F)
ax.plot([R_x0, R_x1], [plot_bottom, plot_bottom],
        transform=fig.transFigure,
        color=NEUTRAL['rule'], linewidth=0.6, zorder=4)
# y-axis caption (rotated, far right)
ax.text(R_x1 + 0.052, (plot_bottom + plot_top) / 2,
        'Count',
        transform=fig.transFigure,
        rotation=270, ha='center', va='center',
        fontsize=9.5, color=NEUTRAL['body'])


# ----- Shared tool legend (top) -----
leg_y = 0.912
total_w = R_x1 - L_x0
swatch_w, swatch_h = 0.014, 0.018
# distribute legend items evenly across the full panel span
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


# ----- Title block -----
ax.text(0.040, 0.972,
        'Fusion-type recall vs false-call profile',
        transform=fig.transFigure,
        ha='left', va='top',
        fontsize=13, fontweight='semibold', color=NEUTRAL['ink'])


save_figure(fig, 'SFig5A', OUTDIR)
plt.close(fig)
