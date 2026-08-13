"""Supplementary Figure 3. Per-tool precision-recall envelopes across the simulation grid.

One facet per tool, ordered by mean F1. Each facet plots the tool's nine
(depth by identity) precision-recall positions, their convex hull as a light
envelope, and the same iso-F1 contours for cross-comparison. The 100 Gb, 95 per
cent identity condition is ringed in every facet as a common anchor.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'SFig3'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

stat = pd.read_csv(DATA / 'stat_summary.tsv', sep='\t')
stat = stat[stat.control == 'positive'].copy()
DEPTH_ORDER = ['1GB', '10GB', '100GB']
ID_ORDER    = ['85%', '90%', '95%']


def collect_points(tool):
    pts, anchor = [], None
    for dep in DEPTH_ORDER:
        for ident in ID_ORDER:
            sub = stat[(stat.Algorithm == tool) & (stat.depth == dep)
                       & (stat.Sequence_Identity == ident)]
            if not len(sub): continue
            p = sub.Precision.iloc[0]; r = sub.Recall.iloc[0]
            if pd.isna(p) or pd.isna(r): continue
            pts.append((float(r), float(p)))
            if dep == '100GB' and ident == '95%':
                anchor = (float(r), float(p))
    return np.array(pts), anchor


tool_data = {tool: collect_points(tool) for tool in TOOL_ORDER}


# ---------- Figure ----------
n = len(TOOL_ORDER)
n_cols = 3
n_rows = 2

fig = plt.figure(figsize=(11.6, 6.4))
fig.patch.set_facecolor('white')

L, R = 0.055, 0.985
B, T = 0.085, 0.770
GAP_X, GAP_Y = 0.030, 0.075
PW = (R - L - GAP_X * (n_cols - 1)) / n_cols
PH = (T - B - GAP_Y * (n_rows - 1)) / n_rows


def _pct(v):
    return '0' if v == 0 else f'{int(round(v*100))}%'


for i, tool in enumerate(TOOL_ORDER):
    row, col = divmod(i, n_cols)
    x = L + col * (PW + GAP_X)
    y = T - (row + 1) * PH - row * GAP_Y
    ax = fig.add_axes([x, y, PW, PH])

    ax.set_xlim(0, 0.32); ax.set_ylim(0, 0.65)
    color = TOOL_COLORS[tool]

    # Iso-F1 contours
    F1_LEVELS = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
    for f1 in F1_LEVELS:
        r_lo = f1 / 2 + 1e-4
        rs = np.linspace(r_lo, 0.40, 2000)
        ps = f1 * rs / (2 * rs - f1)
        mask = (ps >= 0) & (ps <= 0.80) & np.isfinite(ps)
        ax.plot(rs[mask], ps[mask], color='#C5C5C5', linewidth=0.5,
                linestyle=(0, (3, 3)), zorder=1)
    # Labels only in first panel to save space
    if i == 0:
        for f1 in F1_LEVELS:
            LABEL_Y = 0.58
            if 2 * LABEL_Y > f1:
                r_label = LABEL_Y * f1 / (2 * LABEL_Y - f1)
                if 0 < r_label < 0.30:
                    ax.text(r_label + 0.003, LABEL_Y + 0.004,
                            f'{f1:.2f}', fontsize=6.6, color='#9A9A9A',
                            ha='left', va='bottom', zorder=1)

    # Convex hull
    pts, anchor = tool_data[tool]
    if len(pts) >= 3:
        hull = ConvexHull(pts)
        hull_pts = pts[hull.vertices]
        poly = Polygon(hull_pts, facecolor=color, edgecolor=color,
                       alpha=0.30, linewidth=1.4, zorder=3)
        ax.add_patch(poly)
    # Points
    if len(pts):
        ax.scatter(pts[:, 0], pts[:, 1], s=30, c=color,
                   edgecolors='white', linewidths=0.7,
                   alpha=0.95, zorder=4)
    # Anchor — black ring + colored dot
    if anchor is not None:
        ax.scatter(anchor[0], anchor[1], s=120,
                   facecolors='none', edgecolors='black',
                   linewidths=1.8, zorder=5)
        ax.scatter(anchor[0], anchor[1], s=46,
                   c=color, edgecolors='white', linewidths=0.9, zorder=6)

    # Axes styling
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp_name in ('left', 'bottom'):
        ax.spines[sp_name].set_linewidth(1.2)
        ax.spines[sp_name].set_color(NEUTRAL['ink'])
    ax.grid(True, color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
    # Tick label visibility: bottom row gets x ticks, leftmost gets y ticks
    is_left = (col == 0)
    is_bot  = (row == n_rows - 1)
    if is_left:
        ax.set_yticks([0, 0.2, 0.4, 0.6])
        ax.set_yticklabels([_pct(v) for v in [0, 0.2, 0.4, 0.6]],
                           fontsize=8.4)
    else:
        ax.set_yticks([0, 0.2, 0.4, 0.6])
        ax.set_yticklabels([])
    if is_bot:
        ax.set_xticks([0, 0.1, 0.2, 0.3])
        ax.set_xticklabels([_pct(v) for v in [0, 0.1, 0.2, 0.3]],
                           fontsize=8.4)
    else:
        ax.set_xticks([0, 0.1, 0.2, 0.3])
        ax.set_xticklabels([])
    ax.tick_params(axis='both', length=3, width=0.9, color=NEUTRAL['ink'])

    # Panel title with tool color swatch
    ax.set_title(tool, loc='left', fontsize=10.2, fontweight='semibold',
                 color=color, pad=6)


# Shared axis labels
fig.text(L + (R - L) / 2, 0.040, 'Recall',
         ha='center', va='center',
         fontsize=10.6, fontweight='bold', color=NEUTRAL['ink'])
fig.text(0.014, B + (T - B) / 2, 'Precision',
         ha='center', va='center', rotation=90,
         fontsize=10.6, fontweight='bold', color=NEUTRAL['ink'])


# Title
fig.text(0.035, 0.952,
         'Per-tool precision-recall envelopes across the simulation grid',
         ha='left', va='top',
         fontsize=12.4, fontweight='semibold', color=NEUTRAL['ink'])


save_figure(fig, 'SFig3', OUTDIR)
plt.close(fig)
