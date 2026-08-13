"""Supplementary Figure 8B. F1 score in relation to minimum read support.

A three by three grid: rows are sequencing depth (100, 10, 1 Gb) and columns are
read identity (85, 90, 95 per cent). Within each panel the x-axis is the minimum
read-support threshold on a log scale and the y-axis is the F1 score, the
harmonic mean of precision and recall, with one line per caller. Partially
recalled fusions are counted as false positives. This is the per-condition
companion to Fig3A: F1 peaks at an intermediate threshold, and the grid shows how
that optimum shifts with depth and identity.
"""

import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'SFig8'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

METRIC, YLABEL, YLIM = 'F1', 'F1 score', (-0.018, 0.42)
YTICKS = [0, 0.1, 0.2, 0.3, 0.4]
TITLE = 'F1 score in relation to minimum read support'

# ---------- Data ----------
df = pd.read_csv(DATA / 'readsupp_filtering.tsv', sep='\t')
df = df[df.control == 'positive']

DEPTHS = ['100GB', '10GB', '1GB']    # rows
IDENTS = ['85%', '90%', '95%']       # columns
X_MIN, X_MAX = 1, 3000


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
        sub = df[(df.depth == depth) & (df.Sequence_Identity == ident)]

        for tool in reversed(TOOL_ORDER):     # strongest drawn last, on top
            t = sub[sub.Algorithm == tool].sort_values('minimum_read_support')
            if len(t):
                ax.plot(t.minimum_read_support, t[METRIC], color=TOOL_COLORS[tool],
                        linewidth=1.5, alpha=0.92, zorder=3)

        ax.set_xscale('log')
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(*YLIM)
        ax.set_xticks([1, 10, 100, 1000])
        ax.set_yticks(YTICKS)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        for sp in ('left', 'bottom'):
            ax.spines[sp].set_linewidth(1.1)
            ax.spines[sp].set_color(NEUTRAL['ink'])
        ax.grid(True, color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
        ax.tick_params(axis='both', which='both', length=0, labelsize=8.4,
                       color=NEUTRAL['ink'])
        if ci == 0:
            ax.tick_params(axis='y', which='major', length=3.5, width=1.0,
                           color=NEUTRAL['ink'])
        if ri == len(DEPTHS) - 1:
            ax.tick_params(axis='x', which='major', length=3.5, width=1.0,
                           color=NEUTRAL['ink'])

        ax.set_xticklabels(['1', '10', '100', '1000'] if ri == len(DEPTHS) - 1
                           else [])
        ax.set_yticklabels([f'{t:g}' for t in YTICKS] if ci == 0 else [])

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
lab.text((GL + GR) / 2, 0.028, 'Minimum read support', transform=fig.transFigure,
         ha='center', va='center', fontsize=11, fontweight='bold',
         color=NEUTRAL['ink'])
lab.text(0.020, (GB + GT) / 2, YLABEL, transform=fig.transFigure,
         rotation=90, ha='center', va='center', fontsize=11, fontweight='bold',
         color=NEUTRAL['ink'])

# Title
lab.text(0.035, 0.965, TITLE, transform=fig.transFigure, ha='left', va='top',
         fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])

# Tool legend between title and grid, spanning the width.
leg_y = 0.912
leg_left, leg_right = 0.055, 0.945
spacing = (leg_right - leg_left) / len(TOOL_ORDER)
for i, tool in enumerate(TOOL_ORDER):
    lx = leg_left + i * spacing
    lab.scatter(lx, leg_y, s=85, c=TOOL_COLORS[tool], edgecolors='white',
                linewidths=0.7, transform=fig.transFigure, zorder=5, clip_on=False)
    lab.text(lx + 0.011, leg_y, tool, transform=fig.transFigure, ha='left',
             va='center', fontsize=8.5, color=NEUTRAL['body'])


save_figure(fig, 'SFig8B', OUTDIR)
plt.close(fig)
