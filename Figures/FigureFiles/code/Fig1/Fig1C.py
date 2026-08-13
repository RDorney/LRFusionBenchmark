"""Figure 1C. F1 score matrix for the six long-read fusion callers.

For each caller, the F1 score for recovering spiked-in fusions, broken down by
sequencing depth (1Gb, 10Gb, 100Gb) and mean read identity (85%, 90%, 95%).
Rows are ordered by mean F1 at the 100Gb / 95% anchor condition. A right-hand
column gives the per-tool mean across conditions and a bottom row gives the mean
per condition. Cells where a tool produced no output are hatched with an en dash.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig1'
from palette import (
    TOOL_ORDER, TOOL_COLORS, F1_CMAP, NEUTRAL,
    set_style, auto_text_color, save_figure,
)

set_style()
mpl.rcParams['hatch.color'] = '#B5B5B5'
mpl.rcParams['hatch.linewidth'] = 0.45


# ----- Load & pivot -----
df = pd.read_csv(DATA / 'stat_summary.tsv', sep='\t')
df = df[df.control == 'positive'].copy()

DEPTH_ORDER = ['1GB', '10GB', '100GB']
ID_ORDER    = ['85%', '90%', '95%']
COND_ORDER  = [(d, i) for d in DEPTH_ORDER for i in ID_ORDER]

# The source table stores F1 as NA when Precision = Recall = 0 (a divide by zero
# in 2PR/(P+R)). The tool ran and produced zero true calls, so the limiting F1 is
# 0. Set those cells to 0 so they read as zero performance, not missing data.
zero_perf = (df.Precision == 0) & (df.Recall == 0)
df.loc[zero_perf, 'F1'] = 0.0

df['cond'] = df.depth.astype(str) + '|' + df.Sequence_Identity.astype(str)

# Any cell still NaN after the pivot is a condition the tool was never run at.
mat = (
    df.pivot(index='Algorithm', columns='cond', values='F1')
      .reindex(index=TOOL_ORDER,
               columns=[f'{d}|{i}' for d, i in COND_ORDER])
)

vmin, vmax = 0.0, 0.40
row_mean = mat.mean(axis=1, skipna=True)
col_mean = mat.mean(axis=0, skipna=True)

# ----- Layout (in figure-fraction units) -----
fig = plt.figure(figsize=(8.4, 5.0))
fig.patch.set_facecolor('white')

# Heatmap origin and cell size
left    = 0.16
bottom  = 0.20
cell_w  = 0.072
cell_h  = 0.092
n_rows  = len(TOOL_ORDER)
n_cols  = len(COND_ORDER)
gap     = 0.014           # gap between depth groups

# x position of each column
x_pos = []
for j, (d, _i) in enumerate(COND_ORDER):
    extra = DEPTH_ORDER.index(d) * gap
    x_pos.append(left + j * cell_w + extra)
right_edge = x_pos[-1] + cell_w

# Row-mean marginal (right of heatmap)
marg_pad = 0.024
marg_x   = right_edge + marg_pad
marg_w   = cell_w

# Column-mean marginal (below heatmap)
col_marg_pad = 0.018
col_marg_y   = bottom - cell_h - col_marg_pad

# Bracket and label y-positions (above heatmap)
heatmap_top = bottom + n_rows * cell_h
ident_y     = heatmap_top + 0.010          # identity row centre (above heatmap)
bracket_y0  = ident_y + 0.034              # depth bracket bar (top of bracket)
depth_lbl_y = bracket_y0 + 0.013           # depth text just above bracket

# Background / overlay axis (full figure)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_axis_off()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

norm = Normalize(vmin=vmin, vmax=vmax)


# ----- Cell painter -----
def draw_value_cell(x, y, w, h, value, *, font_size=8.6, bold=False):
    """Filled cell with annotated value, auto-contrast text."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        rect = Rectangle((x, y), w, h, transform=fig.transFigure,
                         facecolor=NEUTRAL['na_fill'],
                         edgecolor='white', linewidth=0.9,
                         hatch='////', zorder=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, '—',
                transform=fig.transFigure,
                ha='center', va='center',
                fontsize=11, color=NEUTRAL['na_text'], zorder=3)
        return
    rect = Rectangle((x, y), w, h, transform=fig.transFigure,
                     facecolor=F1_CMAP(norm(value)),
                     edgecolor='white', linewidth=0.9, zorder=2)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, f'{value:.2f}',
            transform=fig.transFigure,
            ha='center', va='center',
            fontsize=font_size,
            color=auto_text_color(value, vmin, vmax, threshold=0.50),
            fontweight=('semibold' if bold else 'normal'),
            zorder=3)


# ----- Main heatmap -----
for r, tool in enumerate(TOOL_ORDER):
    y = bottom + (n_rows - 1 - r) * cell_h
    for c, _ in enumerate(COND_ORDER):
        draw_value_cell(x_pos[c], y, cell_w, cell_h, mat.iloc[r, c])


# ----- Row labels -----
for r, tool in enumerate(TOOL_ORDER):
    y = bottom + (n_rows - 1 - r) * cell_h + cell_h / 2
    ax.text(left - 0.012, y, tool,
            transform=fig.transFigure,
            ha='right', va='center',
            fontsize=9.6, color=NEUTRAL['ink'])


# ----- Column header tier 1: identity (immediately above cells) -----
for c, (_d, ident) in enumerate(COND_ORDER):
    ax.text(x_pos[c] + cell_w/2, ident_y + 0.004, ident,
            transform=fig.transFigure,
            ha='center', va='bottom',
            fontsize=9.0, color=NEUTRAL['body'])


# ----- Column header tier 2: depth group brackets -----
for d in DEPTH_ORDER:
    cols = [j for j, (dd, _i) in enumerate(COND_ORDER) if dd == d]
    x0 = x_pos[cols[0]]
    x1 = x_pos[cols[-1]] + cell_w
    # bracket: horizontal bar with downward tick ends (looks like ┍━━━━━━┑)
    bar_y = bracket_y0
    tick_drop = 0.010
    ax.plot([x0, x1], [bar_y, bar_y], transform=fig.transFigure,
            color=NEUTRAL['rule'], linewidth=0.9,
            solid_capstyle='butt', zorder=3)
    ax.plot([x0, x0], [bar_y, bar_y - tick_drop],
            transform=fig.transFigure,
            color=NEUTRAL['rule'], linewidth=0.9, zorder=3)
    ax.plot([x1, x1], [bar_y, bar_y - tick_drop],
            transform=fig.transFigure,
            color=NEUTRAL['rule'], linewidth=0.9, zorder=3)
    # depth label
    ax.text((x0 + x1) / 2, depth_lbl_y, d,
            transform=fig.transFigure,
            ha='center', va='bottom',
            fontsize=10.2, fontweight='semibold', color=NEUTRAL['ink'])


# ----- Row marginal: per-tool mean across conditions -----
ax.text(marg_x + marg_w/2, ident_y + 0.004, 'mean',
        transform=fig.transFigure,
        ha='center', va='bottom', fontsize=8.6,
        color=NEUTRAL['muted'], style='italic')
for r, tool in enumerate(TOOL_ORDER):
    y = bottom + (n_rows - 1 - r) * cell_h
    v = row_mean.loc[tool]
    rect = Rectangle((marg_x, y), marg_w, cell_h, transform=fig.transFigure,
                     facecolor=F1_CMAP(norm(v)),
                     edgecolor='white', linewidth=0.9, zorder=2)
    ax.add_patch(rect)
    ax.text(marg_x + marg_w/2, y + cell_h/2,
            f'{v:.2f}',
            transform=fig.transFigure,
            ha='center', va='center', fontsize=8.6,
            fontweight='semibold',
            color=auto_text_color(v, vmin, vmax, threshold=0.50), zorder=3)


# ----- Column marginal: per-condition mean across tools -----
ax.text(left - 0.008, col_marg_y + cell_h/2, 'mean',
        transform=fig.transFigure,
        ha='right', va='center', fontsize=8.6,
        color=NEUTRAL['muted'], style='italic')
for c, _ in enumerate(COND_ORDER):
    v = col_mean.iloc[c]
    rect = Rectangle((x_pos[c], col_marg_y), cell_w, cell_h,
                     transform=fig.transFigure,
                     facecolor=F1_CMAP(norm(v)),
                     edgecolor='white', linewidth=0.9, zorder=2)
    ax.add_patch(rect)
    ax.text(x_pos[c] + cell_w/2, col_marg_y + cell_h/2,
            f'{v:.2f}',
            transform=fig.transFigure,
            ha='center', va='center', fontsize=8.6,
            fontweight='semibold',
            color=auto_text_color(v, vmin, vmax, threshold=0.50), zorder=3)


# ----- Colorbar (right of row marginal) -----
cbar_x = marg_x + marg_w + 0.028
cbar_w = 0.013
cbar_h = cell_h * n_rows
ax_cbar = fig.add_axes([cbar_x, bottom, cbar_w, cbar_h])
gradient = np.linspace(0, 1, 256).reshape(-1, 1)
ax_cbar.imshow(gradient, aspect='auto', cmap=F1_CMAP,
               origin='lower', extent=[0, 1, vmin, vmax])
ax_cbar.set_xticks([])
ax_cbar.yaxis.tick_right()
ax_cbar.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
ax_cbar.tick_params(axis='y', length=2.5, pad=2, labelsize=8.4,
                    color=NEUTRAL['rule'])
for spine in ax_cbar.spines.values():
    spine.set_visible(False)
ax_cbar.spines['right'].set_visible(True)
ax_cbar.spines['right'].set_linewidth(0.5)
ax_cbar.spines['right'].set_color(NEUTRAL['rule'])
ax.text(cbar_x + cbar_w + 0.026, bottom + cbar_h / 2, 'F1 score',
        transform=fig.transFigure,
        rotation=90, ha='left', va='center',
        fontsize=9.2, color=NEUTRAL['body'])


# ----- Title -----
ax.text(0.045, 0.955,
        'F1 score across sequencing depth and read identity',
        transform=fig.transFigure,
        ha='left', va='top',
        fontsize=13, fontweight='semibold', color=NEUTRAL['ink'])
ax.text(0.965, 0.955,
        'n = 400 spiked fusions per condition',
        transform=fig.transFigure,
        ha='right', va='top',
        fontsize=8.6, color=NEUTRAL['muted'])


save_figure(fig, 'Fig1C', OUTDIR)
plt.close(fig)
