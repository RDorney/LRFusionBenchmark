"""Figure 3A. Recovery of known and putative novel fusions in K562 and MCF7.

Two panels, one per cell line. Each point is one tool-by-library result on
full-depth SGNex data: x is the novel fusion count (log), y is the known fusion
count. Colour encodes the tool and marker shape encodes the library prep
(circle direct-RNA, triangle direct-cDNA, square PCR-cDNA).
"""

import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig3'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

df = pd.read_csv(DATA / 'counts_SGNex_summary.tsv.gz', sep='\t')
df = df[df.Sequencing_Depth == 'Total'].copy()
# 0 → 0.5 so the log x-axis can show "no novel calls"
df.loc[df.novelnumber == 0, 'novelnumber'] = 0.5

CELL_LINES      = ['K562', 'MCF7']
Y_LIMS          = {'K562': (-0.4, 3.6), 'MCF7': (-2.5, 30)}
LIBRARIES       = ['direct-RNA', 'direct-cDNA', 'PCR-cDNA']
LIBRARY_MARKERS = {'direct-RNA': 'o', 'direct-cDNA': '^', 'PCR-cDNA': 's'}


# ---------- Figure ----------
fig = plt.figure(figsize=(11.4, 6.4))
fig.patch.set_facecolor('white')

L, R = 0.060, 0.985
B, T = 0.120, 0.755
GAP  = 0.075
W = (R - L - GAP) / 2

axes = []
for i, cell in enumerate(CELL_LINES):
    ax = fig.add_axes([L + i * (W + GAP), B, W, T - B])
    axes.append((ax, cell))


for ax, cell in axes:
    ax.set_xscale('log')
    ax.set_xlim(0.3, 1e4)
    ax.set_ylim(*Y_LIMS[cell])

    # Spines / grid
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(1.4)
        ax.spines[sp].set_color(NEUTRAL['ink'])
    ax.grid(True, color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
    ax.tick_params(axis='both', labelsize=9.5, color=NEUTRAL['ink'])

    # Plot points
    sub = df[df.Cell_Lines == cell]
    for _, r in sub.iterrows():
        ax.scatter(r['novelnumber'], r['knownnumber'],
                   s=150, c=TOOL_COLORS[r['Algorithm']],
                   marker=LIBRARY_MARKERS[r['Library']],
                   edgecolors='black', linewidths=1.0,
                   alpha=0.88, zorder=4, clip_on=False)

    ax.set_title(cell, fontsize=12.4, fontweight='semibold',
                 color=NEUTRAL['ink'], loc='left', pad=10)
    ax.set_xlabel('Novel fusions  (count, log)',
                  fontsize=10.6, fontweight='bold',
                  color=NEUTRAL['ink'], labelpad=4)
    ax.set_ylabel('Known fusions recovered',
                  fontsize=10.6, fontweight='bold',
                  color=NEUTRAL['ink'], labelpad=4)


# ---------- Title block + legends ----------
ax_lab = fig.add_axes([0, 0, 1, 1])
ax_lab.set_axis_off()
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)

ax_lab.text(0.035, 0.962,
            'Recovery of known fusions and putative novel fusions in '
            'K562 and MCF7 cells',
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])

# Tool colour legend (scatter dots, evenly spaced)
leg_y = 0.890
leg_left, leg_right = 0.040, 0.660
spacing = (leg_right - leg_left) / len(TOOL_ORDER)
for i, tool in enumerate(TOOL_ORDER):
    lx = leg_left + i * spacing + 0.005
    ax_lab.scatter(lx, leg_y, s=90, c=TOOL_COLORS[tool],
                   edgecolors='white', linewidths=0.7,
                   transform=fig.transFigure, zorder=4, clip_on=False)
    ax_lab.text(lx + 0.012, leg_y, tool,
                transform=fig.transFigure,
                ha='left', va='center',
                fontsize=8.7, color=NEUTRAL['body'])

# Library shape legend (right of tool legend)
shape_left = 0.700
shape_spacing = 0.095
for i, lib in enumerate(LIBRARIES):
    sx = shape_left + i * shape_spacing
    ax_lab.scatter(sx, leg_y, s=90, marker=LIBRARY_MARKERS[lib],
                   c='#888', edgecolors='black', linewidths=0.8,
                   transform=fig.transFigure, zorder=4, clip_on=False)
    ax_lab.text(sx + 0.012, leg_y, lib,
                transform=fig.transFigure,
                ha='left', va='center',
                fontsize=8.7, color=NEUTRAL['body'])


save_figure(fig, 'Fig3A', OUTDIR)
plt.close(fig)
