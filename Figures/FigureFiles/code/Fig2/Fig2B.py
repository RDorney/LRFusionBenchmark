"""Figure 2B. Effect of the minimum read-support threshold on F1, precision and recall.

Three stacked panels, one per metric, with the six callers overlaid as lines
against the minimum read-support threshold on a log x-axis. Conditions: 100Gb
depth, 95% mean identity, positive control.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig2'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

df = pd.read_csv(DATA / 'readsupp_filtering.tsv', sep='\t')
df = df[(df.depth == '100GB') & (df.Sequence_Identity == '95%')
        & (df.control == 'positive')].copy()
df = df.sort_values(['Algorithm', 'minimum_read_support'])


# ---------- Figure ----------
fig = plt.figure(figsize=(8.4, 8.0))
fig.patch.set_facecolor('white')

# Three stacked panels, sharing x-axis
LEFT, RIGHT = 0.115, 0.965
PANEL_BOTTOM = 0.090
PANEL_TOP    = 0.815
TOTAL_H = PANEL_TOP - PANEL_BOTTOM
GAP = 0.020
PANEL_H = (TOTAL_H - 2 * GAP) / 3

ax_F1 = fig.add_axes([LEFT, PANEL_TOP - PANEL_H,                       RIGHT - LEFT, PANEL_H])
ax_P  = fig.add_axes([LEFT, PANEL_TOP - 2 * PANEL_H - GAP,             RIGHT - LEFT, PANEL_H])
ax_R  = fig.add_axes([LEFT, PANEL_TOP - 3 * PANEL_H - 2 * GAP,         RIGHT - LEFT, PANEL_H])

X_MIN, X_MAX = 1, 3000


def style_panel(ax, ylabel, ylim, show_x=False):
    ax.set_xscale('log')
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontsize=10.5, fontweight='bold', color=NEUTRAL['ink'])
    ax.grid(True, color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
    for spine in ax.spines.values():
        spine.set_color(NEUTRAL['rule'])
        spine.set_linewidth(0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.4)
    ax.spines['left'].set_color(NEUTRAL['ink'])
    if show_x:
        ax.spines['bottom'].set_linewidth(1.4)
        ax.spines['bottom'].set_color(NEUTRAL['ink'])
        ax.tick_params(axis='x', labelsize=9, length=4, width=1.0, color=NEUTRAL['ink'])
    else:
        ax.tick_params(axis='x', labelsize=0, length=0)
    ax.tick_params(axis='y', labelsize=9, length=4, width=1.0, color=NEUTRAL['ink'])


style_panel(ax_F1, 'F1',        (-0.018, 0.42))
style_panel(ax_P,  'Precision', (-0.044, 1.05))
style_panel(ax_R,  'Recall',    (-0.014, 0.32), show_x=True)
ax_R.set_xlabel('Minimum read support', fontsize=11, fontweight='bold',
                color=NEUTRAL['ink'], labelpad=4)


# ---------- Draw 6 tool lines per panel ----------
draw_order = list(reversed(TOOL_ORDER))   # weakest first so strongest on top

for tool in draw_order:
    t = df[df.Algorithm == tool].sort_values('minimum_read_support')
    if not len(t):
        continue
    color = TOOL_COLORS[tool]
    rs = t.minimum_read_support.values
    f1 = t.F1.values
    p  = t.Precision.values
    r  = t.Recall.values

    ax_F1.plot(rs, f1, color=color, linewidth=1.8, zorder=3, alpha=0.92)
    ax_P.plot(rs, p,  color=color, linewidth=1.8, zorder=3, alpha=0.92)
    ax_R.plot(rs, r,  color=color, linewidth=1.8, zorder=3, alpha=0.92)


# ---------- Title + legend ----------
ax_lab = fig.add_axes([0, 0, 1, 1])
ax_lab.set_axis_off()
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)

ax_lab.text(0.035, 0.965,
            'Effect of minimum read-support threshold on F1, Precision, Recall',
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=12.5, fontweight='semibold', color=NEUTRAL['ink'])
# Tool legend
leg_y = 0.880
leg_left, leg_right = 0.040, 0.965
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


save_figure(fig, 'Fig2B', OUTDIR)
plt.close(fig)
