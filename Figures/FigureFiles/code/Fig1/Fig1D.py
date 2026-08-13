"""Figure 1D. Precision-recall position across the simulation grid.

A three by three grid of precision-recall scatter plots. Rows are mean read
identity (95, 90 and 85 per cent, top to bottom) and columns are sequencing
depth (1, 10 and 100 Gb, left to right). Each cell shows the six tools as dots
over a background of iso-F1 contours, with the highest-F1 tool at that condition
drawn with a thicker ink ring. Only the leftmost column and bottom row carry
tick numbers, keeping the inner panels to contours and data.
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


# ---------- Data ----------
stat = pd.read_csv(DATA / 'stat_summary.tsv', sep='\t')
stat = stat[stat.control == 'positive'].copy()


def get(metric, tool, depth, identity):
    sub = stat[(stat.Algorithm == tool) & (stat.depth == depth)
               & (stat.Sequence_Identity == identity)]
    if not len(sub): return None
    v = sub[metric].iloc[0]
    return None if pd.isna(v) else float(v)


DEPTH_COLS = ['1GB', '10GB', '100GB']
ID_ROWS    = ['95%', '90%', '85%']

X_LIM = 0.32
Y_LIM = 0.65
F1_LEVELS = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50]


def best_tool_at(depth, identity):
    """Return the tool with highest F1 at this condition (None if all NaN)."""
    best_tool, best_f1 = None, -1.0
    for tool in TOOL_ORDER:
        p = get('Precision', tool, depth, identity)
        r = get('Recall',    tool, depth, identity)
        if p is None or r is None: continue
        if (p + r) == 0:
            f1 = 0.0
        else:
            f1 = 2 * p * r / (p + r)
        if f1 > best_f1:
            best_f1 = f1
            best_tool = tool
    return best_tool


# ---------- Figure ----------
fig = plt.figure(figsize=(11.0, 8.6))
fig.patch.set_facecolor('white')

GRID_left, GRID_right = 0.078, 0.795
GRID_bottom, GRID_top = 0.090, 0.785
gap_x, gap_y = 0.012, 0.014
n_rows = n_cols = 3
panel_w = (GRID_right - GRID_left - (n_cols - 1) * gap_x) / n_cols
panel_h = (GRID_top    - GRID_bottom - (n_rows - 1) * gap_y) / n_rows


def panel_axes(r_idx, c_idx):
    x0 = GRID_left + c_idx * (panel_w + gap_x)
    y0 = GRID_top  - (r_idx + 1) * panel_h - r_idx * gap_y
    return fig.add_axes([x0, y0, panel_w, panel_h])


# Overlay for global labels
ax_lab = fig.add_axes([0, 0, 1, 1])
ax_lab.set_axis_off()
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)


# ---------- Render each panel ----------
def short_tick(v):
    """Format a 0–1 fraction as a short percentage: 0 → '0', 0.1 → '10%'."""
    if v == 0:
        return '0'
    return f'{int(round(v * 100))}%'


for r_idx, ident in enumerate(ID_ROWS):
    for c_idx, dep in enumerate(DEPTH_COLS):
        ax = panel_axes(r_idx, c_idx)
        ax.set_xlim(0, X_LIM)
        ax.set_ylim(0, Y_LIM)

        is_leftmost = (c_idx == 0)
        is_bottom   = (r_idx == n_rows - 1)

        # iso-F1 contours with exact-position numeric labels on every panel.
        for f1 in F1_LEVELS:
            r_lo = f1 / 2 + 1e-4
            rs = np.linspace(r_lo, 0.40, 2000)
            ps = f1 * rs / (2 * rs - f1)
            mask = (ps >= 0) & (ps <= 0.80) & np.isfinite(ps)
            ax.plot(rs[mask], ps[mask],
                    color='#9C9C9C', linewidth=1.0,
                    linestyle=(0, (4, 3)),
                    zorder=1)
            # Exact-position labels at LABEL_Y near the top of each panel
            LABEL_Y_GRID = Y_LIM - 0.07
            if 2 * LABEL_Y_GRID > f1:
                r_label = LABEL_Y_GRID * f1 / (2 * LABEL_Y_GRID - f1)
                if 0 < r_label < X_LIM - 0.02:
                    ax.text(r_label + 0.003, LABEL_Y_GRID + 0.005,
                            f'{f1:.2f}',
                            fontsize=6.4, color='#7A7A7A',
                            ha='left', va='bottom', zorder=2)

        # Determine the best tool at this condition for the highlight ring
        best = best_tool_at(dep, ident)

        # Plot dots
        for tool in TOOL_ORDER:
            p = get('Precision', tool, dep, ident)
            r = get('Recall',    tool, dep, ident)
            if p is None or r is None:
                continue
            is_best = (tool == best)
            color = TOOL_COLORS[tool]
            if is_best:
                # ink-colour outer ring, then the colour-filled inner dot on top
                ax.scatter(r, p, s=220, c=NEUTRAL['ink'],
                           zorder=4, clip_on=False)
                ax.scatter(r, p, s=130, c=color,
                           edgecolors='white', linewidths=0.9,
                           zorder=5, clip_on=False)
            else:
                ax.scatter(r, p, s=130, c=color,
                           edgecolors='white', linewidths=1.0,
                           zorder=4, clip_on=False)

        # (No per-panel inset condition label — the outer row & column headers
        #  already say which condition each cell represents.)

        # Spine + tick visibility
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

        if is_leftmost and is_bottom:
            # both axes shown, bold
            for side in ('left', 'bottom'):
                ax.spines[side].set_visible(True)
                ax.spines[side].set_linewidth(1.6)
                ax.spines[side].set_color(NEUTRAL['ink'])
            ax.tick_params(axis='x', labelsize=8.4, length=4, width=1.2,
                           color=NEUTRAL['ink'])
            ax.tick_params(axis='y', labelsize=8.4, length=4, width=1.2,
                           color=NEUTRAL['ink'])
            # short tick labels
            xt = [0, 0.1, 0.2, 0.3]
            yt = [0, 0.2, 0.4, 0.6]
            ax.set_xticks(xt); ax.set_xticklabels([short_tick(t) for t in xt])
            ax.set_yticks(yt); ax.set_yticklabels([short_tick(t) for t in yt])
        elif is_leftmost:
            ax.spines['left'].set_visible(True)
            ax.spines['left'].set_linewidth(1.6)
            ax.spines['left'].set_color(NEUTRAL['ink'])
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(axis='y', labelsize=8.4, length=4, width=1.2,
                           color=NEUTRAL['ink'])
            ax.tick_params(axis='x', length=0, labelsize=0)
            yt = [0, 0.2, 0.4, 0.6]
            ax.set_yticks(yt); ax.set_yticklabels([short_tick(t) for t in yt])
            ax.set_xticks([])
        elif is_bottom:
            ax.spines['bottom'].set_visible(True)
            ax.spines['bottom'].set_linewidth(1.6)
            ax.spines['bottom'].set_color(NEUTRAL['ink'])
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis='x', labelsize=8.4, length=4, width=1.2,
                           color=NEUTRAL['ink'])
            ax.tick_params(axis='y', length=0, labelsize=0)
            xt = [0, 0.1, 0.2, 0.3]
            ax.set_xticks(xt); ax.set_xticklabels([short_tick(t) for t in xt])
            ax.set_yticks([])
        else:
            # inner panels: no spines, no ticks, no gridlines
            for side in ('left', 'bottom'):
                ax.spines[side].set_visible(False)
            ax.set_xticks([]); ax.set_yticks([])
        # never any background gridlines
        ax.grid(False)


# ---------- Outer axis labels (bold) ----------
ax_lab.text(GRID_left - 0.045, (GRID_bottom + GRID_top) / 2,
            'Precision',
            transform=fig.transFigure,
            rotation=90, ha='center', va='center',
            fontsize=11.5, fontweight='bold', color=NEUTRAL['ink'])
ax_lab.text((GRID_left + GRID_right) / 2, GRID_bottom - 0.045,
            'Recall',
            transform=fig.transFigure,
            ha='center', va='top',
            fontsize=11.5, fontweight='bold', color=NEUTRAL['ink'])


# ---------- Outer column headers (depth) ----------
header_y = GRID_top + 0.018
for c_idx, dep in enumerate(DEPTH_COLS):
    cx = GRID_left + (c_idx + 0.5) * panel_w + c_idx * gap_x
    ax_lab.text(cx, header_y, dep,
                transform=fig.transFigure,
                ha='center', va='bottom',
                fontsize=11.5, fontweight='semibold', color=NEUTRAL['ink'])
ax_lab.text((GRID_left + GRID_right) / 2, header_y + 0.028,
            'Sequencing depth',
            transform=fig.transFigure,
            ha='center', va='bottom',
            fontsize=9.5, color=NEUTRAL['muted'], style='italic')


# ---------- Outer row headers (identity) on right side ----------
for r_idx, ident in enumerate(ID_ROWS):
    cy = GRID_top - (r_idx + 0.5) * panel_h - r_idx * gap_y
    ax_lab.text(GRID_right + 0.012, cy, ident,
                transform=fig.transFigure,
                ha='left', va='center',
                fontsize=11.5, fontweight='semibold', color=NEUTRAL['ink'])
ax_lab.text(GRID_right + 0.050, (GRID_bottom + GRID_top) / 2,
            'Read identity',
            transform=fig.transFigure,
            rotation=270, ha='center', va='center',
            fontsize=9.5, color=NEUTRAL['muted'], style='italic')


# ---------- Tool legend (top row) — lowered a bit to clear the subtitle ----------
leg_y = 0.880
total_w = GRID_right - GRID_left
spacing = total_w / len(TOOL_ORDER)
for i, tool in enumerate(TOOL_ORDER):
    lx = GRID_left + i * spacing + 0.008
    ax_lab.scatter(lx, leg_y, s=110, c=TOOL_COLORS[tool],
                   edgecolors='white', linewidths=0.8,
                   transform=fig.transFigure, zorder=4, clip_on=False)
    ax_lab.text(lx + 0.012, leg_y, tool,
                transform=fig.transFigure,
                ha='left', va='center',
                fontsize=8.8, color=NEUTRAL['body'])

# (F1 numeric labels live inside the top-left panel directly on the contours.)


# ---------- Title ----------
ax_lab.text(0.025, 0.965,
            'Precision-recall position across the simulation grid',
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=13, fontweight='semibold', color=NEUTRAL['ink'])


save_figure(fig, 'Fig1D', OUTDIR)
plt.close(fig)
