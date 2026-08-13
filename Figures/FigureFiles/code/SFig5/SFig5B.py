"""Supplementary Figure 5B. Partially recalled fusions at 100 Gb and 95% identity.

Panel A gives the total number of partially recalled fusions per tool. Panel B
breaks that total into the eight miss-call subtypes, one facet each, with a bar
per tool. A partial recall is a spiked fusion that was detected but with the
wrong structure: CM = chromosomal misalignment, R = reversed gene order, W =
wrong order (tri-fusions only), and truncated tri-fusions where a gene is lost.
Bars are coloured by tool. The condition is fixed at 100 Gb depth and 95% read
identity, where recall is highest.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'SFig5'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

# ---------- Data (100 Gb, 95% identity, positive control, partial recalls) ----------
df = pd.read_csv(DATA / 'fusion_profile_stat_summary.tsv', sep='\t')
df = df[(df.control == 'positive') & (df.depth == '100GB')
        & (df.Sequence_Identity == '95%')
        & (df.recall_category == 'Partial_Recall')]

total_per_tool = df.groupby('Algorithm')['fusion_type_count'].sum()
by_subtype = (df.groupby(['fusionType', 'Algorithm'])['fusion_type_count'].sum()
              .unstack(fill_value=0).reindex(columns=TOOL_ORDER, fill_value=0))

# Subtypes ordered by total, with display labels.
SUBTYPES = [
    ('truncated_tri_fusion',                       'Truncated\nTri-fusion'),
    ('reverse_order:truncated_tri_fusion',         '(R) Truncated\nTri-fusion'),
    ('reverse_order:intra_chromosome',             '(R) Intra-\nchromosomal'),
    ('reverse_order:hybrid',                       '(R) Inter-\nchromosomal'),
    ('reverse_order:read_through',                 '(R) Read-through'),
    ('wrong_order:tri_fusion',                     '(W) Tri-fusion'),
    ('chromosomal_misalignment:intra_chromosome',  '(CM) Intra-\nchromosomal'),
    ('reverse_order:tri_fusion',                   '(R) Tri-fusion'),
]

bar_colors = [TOOL_COLORS[t] for t in TOOL_ORDER]


def style_axes(ax):
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(1.1)
        ax.spines[sp].set_color(NEUTRAL['ink'])
    ax.grid(True, axis='y', color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
    ax.tick_params(axis='both', length=0, labelsize=8.4, color=NEUTRAL['ink'])


# ---------- Figure ----------
fig = plt.figure(figsize=(13.6, 6.2))
fig.patch.set_facecolor('white')

# Panel A: total per tool
axA = fig.add_axes([0.058, 0.135, 0.185, 0.655])
axA.bar(range(len(TOOL_ORDER)), [total_per_tool.get(t, 0) for t in TOOL_ORDER],
        color=bar_colors, edgecolor='white', linewidth=0.7, zorder=3)
axA.set_xticks(range(len(TOOL_ORDER)))
axA.set_xticklabels(TOOL_ORDER, rotation=40, ha='right', fontsize=8.4,
                    color=NEUTRAL['body'])
axA.set_ylim(0, 200)
axA.set_ylabel('Partially recalled fusions', fontsize=10.4, fontweight='bold',
               color=NEUTRAL['ink'])
style_axes(axA)

# Panel B: 2 x 4 facet grid, one per subtype
B_left, B_right = 0.340, 0.988
B_bottom, B_top = 0.135, 0.790
gap_x, gap_y = 0.020, 0.150
bw = (B_right - B_left - 3 * gap_x) / 4
bh = (B_top - B_bottom - gap_y) / 2

for i, (ft, label) in enumerate(SUBTYPES):
    row, col = divmod(i, 4)
    x0 = B_left + col * (bw + gap_x)
    y0 = B_top - (row + 1) * bh - row * gap_y
    ax = fig.add_axes([x0, y0, bw, bh])
    vals = [by_subtype.loc[ft, t] if ft in by_subtype.index else 0
            for t in TOOL_ORDER]
    ax.bar(range(len(TOOL_ORDER)), vals, color=bar_colors,
           edgecolor='white', linewidth=0.6, zorder=3)
    ax.set_ylim(0, 100)
    ax.set_xticks([])
    style_axes(ax)
    if col == 0:
        ax.set_yticks([0, 50, 100])
    else:
        ax.set_yticks([0, 50, 100])
        ax.set_yticklabels([])
    ax.set_title(label, fontsize=8.6, color=NEUTRAL['ink'], pad=4,
                 linespacing=0.95)


# ---------- Overlay: title, tags, legend, note ----------
lab = fig.add_axes([0, 0, 1, 1])
lab.set_axis_off(); lab.set_xlim(0, 1); lab.set_ylim(0, 1)

lab.text(0.035, 0.965, 'Partially recalled fusions  (100 Gb, 95% identity)',
         transform=fig.transFigure, ha='left', va='top',
         fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])

lab.text(0.058, 0.876, 'Total per tool', transform=fig.transFigure,
         ha='left', va='top', fontsize=9.5, color=NEUTRAL['muted'])
lab.text(0.340, 0.876, 'By miss-call subtype', transform=fig.transFigure,
         ha='left', va='top', fontsize=9.5, color=NEUTRAL['muted'])

# Tool legend (top right)
leg_x = 0.560
for tool in TOOL_ORDER:
    lab.scatter(leg_x, 0.955, s=80, c=TOOL_COLORS[tool], edgecolors='white',
                linewidths=0.7, transform=fig.transFigure, zorder=5,
                clip_on=False)
    lab.text(leg_x + 0.008, 0.955, tool, transform=fig.transFigure,
             ha='left', va='center', fontsize=8.2, color=NEUTRAL['body'])
    leg_x += 0.011 + 0.0088 * len(tool)

# Category key
lab.text(0.340, 0.045,
         'CM = chromosomal misalignment    R = reversed order    '
         'W = wrong order (tri-fusions only)',
         transform=fig.transFigure, ha='left', va='top',
         fontsize=8.2, color=NEUTRAL['muted'])


save_figure(fig, 'SFig5B', OUTDIR)
plt.close(fig)
