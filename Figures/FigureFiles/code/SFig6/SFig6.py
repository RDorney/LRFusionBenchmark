"""Supplementary Figure 6. Breakpoint-position accuracy across the simulation grid.

A three by three grid: rows are sequencing depth (100, 10, 1 Gb) and columns are
read identity (85, 90, 95 per cent). Within each panel the callers are grouped
into Breakpoint 1 and Breakpoint 2, with one boxplot per caller of the absolute
distance between its predicted breakpoint and the simulated truth on a log
y-axis. The dashed line marks 100 bp. Low-depth panels are sparse because few
fusions are recalled there.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'SFig6'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

# ---------- Data ----------
df = pd.read_csv(DATA / 'absolute_Breakpoint_Accuracy.tsv', sep='\t')
df = df[df.Breakpoint_number.isin(['Breakpoint_1', 'Breakpoint_2'])].copy()
df['abs_d'] = pd.to_numeric(df['Distance_from_Simulated_Breakpoint'],
                            errors='coerce').abs()
df = df.dropna(subset=['abs_d'])
df.loc[df.abs_d == 0, 'abs_d'] = 0.5    # 0 bp plotted at the log floor

DEPTHS = ['100GB', '10GB', '1GB']    # rows
IDENTS = ['85%', '90%', '95%']       # columns
Y_MIN, Y_MAX = 0.4, 5e8
REF_100 = 100
YTICKS = [1, 100, 1e4, 1e6, 1e8]
YLABELS = ['1', '$10^2$', '$10^4$', '$10^6$', '$10^8$']


# ---------- Figure ----------
fig = plt.figure(figsize=(10.6, 8.6))
fig.patch.set_facecolor('white')

GL, GR = 0.088, 0.912
GT, GB = 0.800, 0.095
gx, gy = 0.028, 0.052
pw = (GR - GL - 2 * gx) / 3
ph = (GT - GB - 2 * gy) / 3

rng = np.random.default_rng(0)

BPS = [('Breakpoint_1', 'Breakpoint 1'), ('Breakpoint_2', 'Breakpoint 2')]
N_TOOLS = len(TOOL_ORDER)
GROUP = N_TOOLS + 1        # a blank slot separates the two breakpoint blocks

for ri, depth in enumerate(DEPTHS):
    for ci, ident in enumerate(IDENTS):
        x0 = GL + ci * (pw + gx)
        y0 = GT - (ri + 1) * ph - ri * gy
        ax = fig.add_axes([x0, y0, pw, ph])
        sub = df[(df.seq_depth == depth) & (df.seq_id == ident)]

        ax.set_yscale('log')
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.set_xlim(0.3, 2 * GROUP - 0.3)
        ax.axhline(REF_100, color='#B0B7BE', linewidth=1.0,
                   linestyle=(0, (3, 2)), zorder=1)

        box_data, positions, box_colors = [], [], []
        for bi, (bp_name, _lbl) in enumerate(BPS):
            for ti, tool in enumerate(TOOL_ORDER):
                v = sub[(sub.Breakpoint_number == bp_name)
                        & (sub.Algorithm == tool)]['abs_d'].values
                pos = bi * GROUP + ti + 1
                box_data.append(v if len(v) else [np.nan])
                positions.append(pos)
                box_colors.append(TOOL_COLORS[tool])
                if len(v):
                    ax.scatter(pos + rng.uniform(-0.2, 0.2, len(v)), v, s=4,
                               c=[TOOL_COLORS[tool]], edgecolors=NEUTRAL['ink'],
                               linewidths=0.2, alpha=0.55, zorder=4)

        bx = ax.boxplot(box_data, positions=positions, widths=0.72,
                        whis=(0, 100), showfliers=False, patch_artist=True,
                        zorder=3,
                        medianprops=dict(color=NEUTRAL['ink'], linewidth=1.4),
                        whiskerprops=dict(color=NEUTRAL['ink'], linewidth=0.8),
                        capprops=dict(color=NEUTRAL['ink'], linewidth=0.8))
        for i, c in enumerate(box_colors):
            bx['boxes'][i].set(facecolor=c, alpha=0.85, edgecolor=NEUTRAL['ink'],
                               linewidth=0.8)

        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        for sp in ('left', 'bottom'):
            ax.spines[sp].set_linewidth(1.1)
            ax.spines[sp].set_color(NEUTRAL['ink'])
        ax.grid(True, axis='y', color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
        ax.set_xticks([])
        ax.set_yticks(YTICKS)
        ax.tick_params(axis='both', which='both', length=0, labelsize=8.4,
                       color=NEUTRAL['ink'])
        if ci == 0:
            ax.set_yticklabels(YLABELS)
            ax.tick_params(axis='y', which='major', length=3.5, width=1.0,
                           color=NEUTRAL['ink'])
        else:
            ax.set_yticklabels([])

        # Breakpoint block labels under the bottom row
        if ri == len(DEPTHS) - 1:
            for bi, (_bp, lbl) in enumerate(BPS):
                ax.text(bi * GROUP + (N_TOOLS + 1) / 2, -0.035, lbl,
                        transform=ax.get_xaxis_transform(), ha='center',
                        va='top', fontsize=8.8, color=NEUTRAL['body'])

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
lab.text(0.022, (GB + GT) / 2,
         '|Distance from simulated breakpoint|  (bp, log)',
         transform=fig.transFigure, rotation=90, ha='center', va='center',
         fontsize=10.6, fontweight='bold', color=NEUTRAL['ink'])

# Title
lab.text(0.035, 0.965, 'Breakpoint-position accuracy across the simulation grid',
         transform=fig.transFigure, ha='left', va='top',
         fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])

# Tool legend between title and grid
leg_y = 0.912
leg_left, leg_right = 0.055, 0.900
spacing = (leg_right - leg_left) / len(TOOL_ORDER)
for i, tool in enumerate(TOOL_ORDER):
    lx = leg_left + i * spacing
    lab.scatter(lx, leg_y, s=85, c=TOOL_COLORS[tool], edgecolors='white',
                linewidths=0.7, transform=fig.transFigure, zorder=5, clip_on=False)
    lab.text(lx + 0.011, leg_y, tool, transform=fig.transFigure, ha='left',
             va='center', fontsize=8.5, color=NEUTRAL['body'])
lab.text(0.906, leg_y, '- - 100 bp', transform=fig.transFigure, ha='left',
         va='center', fontsize=8.2, color=NEUTRAL['muted'])


save_figure(fig, 'SFig6', OUTDIR)
plt.close(fig)
