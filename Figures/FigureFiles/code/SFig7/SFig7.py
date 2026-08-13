"""Supplementary Figure 7. Unique fusions by number of supporting algorithms.

A three by three grid of line plots: rows are sequencing depth (100, 10, 1 Gb)
and columns are read identity (85, 90, 95 per cent). Within each panel the x-axis
is the number of algorithms that called a fusion (1 to 6) and the log y-axis is
the number of unique fusions, split into three lines by recall category (false
calls, partial recalls, true recalls). Colours match SFig5. The pattern is
consistent across the grid: false calls concentrate at a single algorithm and
fall away quickly, while true fusions persist across many algorithms, and the
depth of agreement grows with sequencing depth.
"""

import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'SFig7'
from palette import PRECISION_CMAP, RECALL_CMAP, F1_CMAP, NEUTRAL, set_style, save_figure

set_style()

# ---------- Data ----------
df = pd.read_csv(DATA / 'fusions_in_n_algorithms.tsv', sep='\t')
df = df[df.control == 'positive']
uniq = df.drop_duplicates(['fusion.gene.id', 'overlap', 'recall_category',
                           'depth', 'Sequence_Identity'])

DEPTHS = ['100GB', '10GB', '1GB']    # rows, top to bottom
IDENTS = ['85%', '90%', '95%']       # columns, left to right
CATS = [
    ('False_Call',     'False call',     PRECISION_CMAP(0.78)),
    ('Partial_Recall', 'Partial recall', RECALL_CMAP(0.72)),
    ('True_Recall',    'True recall',    F1_CMAP(0.62)),
]


def summary(depth, ident):
    sub = uniq[(uniq.depth == depth) & (uniq.Sequence_Identity == ident)]
    return (sub.groupby(['overlap', 'recall_category']).size()
            .unstack(fill_value=0).reindex(range(1, 7), fill_value=0)
            .reindex(columns=[c for c, _, _ in CATS], fill_value=0))


# ---------- Figure ----------
fig = plt.figure(figsize=(10.6, 8.6))
fig.patch.set_facecolor('white')

GL, GR = 0.088, 0.915
GT, GB = 0.800, 0.095
gx, gy = 0.028, 0.052
pw = (GR - GL - 2 * gx) / 3
ph = (GT - GB - 2 * gy) / 3

for ri, depth in enumerate(DEPTHS):
    for ci, ident in enumerate(IDENTS):
        x0 = GL + ci * (pw + gx)
        y0 = GT - (ri + 1) * ph - ri * gy
        ax = fig.add_axes([x0, y0, pw, ph])
        s = summary(depth, ident)

        for col, _label, color in CATS:
            ys = s[col]
            mask = ys > 0
            ax.plot(s.index[mask], ys[mask], marker='o', markersize=5,
                    linewidth=1.8, color=color, markeredgecolor='white',
                    markeredgewidth=0.7, zorder=3, clip_on=False)

        ax.set_yscale('log')
        ax.set_ylim(0.8, 4000)
        ax.set_xlim(0.7, 6.3)
        ax.set_xticks(range(1, 7))
        ax.set_yticks([1, 10, 100, 1000])
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        for sp in ('left', 'bottom'):
            ax.spines[sp].set_linewidth(1.1)
            ax.spines[sp].set_color(NEUTRAL['ink'])
        ax.grid(True, axis='y', color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
        ax.tick_params(axis='both', which='both', length=0, labelsize=8.4,
                       color=NEUTRAL['ink'])
        if ci == 0:
            ax.tick_params(axis='y', which='major', length=3.5, width=1.0,
                           color=NEUTRAL['ink'])
        if ri == len(DEPTHS) - 1:
            ax.tick_params(axis='x', which='major', length=3.5, width=1.0,
                           color=NEUTRAL['ink'])

        if ri == len(DEPTHS) - 1:
            ax.set_xticklabels(range(1, 7))
        else:
            ax.set_xticklabels([])
        if ci == 0:
            ax.set_yticklabels(['1', '10', '100', '1000'])
        else:
            ax.set_yticklabels([])

        if ri == 0:
            ax.set_title(ident, fontsize=11, fontweight='semibold',
                         color=NEUTRAL['ink'], pad=8)
        if ci == len(IDENTS) - 1:
            ax.text(1.045, 0.5, depth, transform=ax.transAxes, rotation=270,
                    ha='left', va='center', fontsize=11, fontweight='semibold',
                    color=NEUTRAL['ink'])


# ---------- Outer labels ----------
lab = fig.add_axes([0, 0, 1, 1])
lab.set_axis_off(); lab.set_xlim(0, 1); lab.set_ylim(0, 1)

lab.text((GL + GR) / 2, GT + 0.052, 'Read identity', transform=fig.transFigure,
         ha='center', va='bottom', fontsize=9.5, color=NEUTRAL['muted'],
         style='italic')
lab.text(GR + 0.052, (GB + GT) / 2, 'Sequencing depth', transform=fig.transFigure,
         rotation=270, ha='center', va='center', fontsize=9.5,
         color=NEUTRAL['muted'], style='italic')

lab.text((GL + GR) / 2, 0.028, 'Supporting algorithms', transform=fig.transFigure,
         ha='center', va='center', fontsize=11, fontweight='bold',
         color=NEUTRAL['ink'])
lab.text(0.020, (GB + GT) / 2, 'Number of unique fusions',
         transform=fig.transFigure, rotation=90, ha='center', va='center',
         fontsize=11, fontweight='bold', color=NEUTRAL['ink'])

# Title + legend
lab.text(0.035, 0.965,
         'Unique fusions by number of supporting algorithms',
         transform=fig.transFigure, ha='left', va='top',
         fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])

# Legend between the title and the grid, centred, matching the tool legends.
LEG_Y = 0.910
widths = [0.028 + 0.0075 * len(label) for _col, label, _color in CATS]
leg_x = (GL + GR) / 2 - sum(widths) / 2
for (_col, label, color), w in zip(CATS, widths):
    lab.scatter(leg_x, LEG_Y, s=90, c=[color], edgecolors='white',
                linewidths=0.7, transform=fig.transFigure, zorder=5,
                clip_on=False)
    lab.text(leg_x + 0.011, LEG_Y, label, transform=fig.transFigure,
             ha='left', va='center', fontsize=9.2, color=NEUTRAL['body'])
    leg_x += w


save_figure(fig, 'SFig7', OUTDIR)
plt.close(fig)
