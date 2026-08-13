"""Supplementary Figure 2. Precision, recall and F1 with partial calls scored as true recall.

Recalculated counterpart to SFig1. A called fusion is classed as true recall,
partial recall (truncated, reversed, chromosomally misaligned or wrong order) or
false call. SFig1 counts partial calls as false; here they are counted as true
recall instead, moving them from the false pool into the true pool before the
metrics are recomputed:

    true positives  = true recall + partial recall
    false positives = false calls
    recall    = true positives / spiked fusions
    precision = true positives / (true positives + false positives)

Layout matches SFig1: three stacked heatmaps (one per metric), six tools ordered
by mean F1, nine conditions (three depths x three identities) under a two-tier
header, each metric on its own colour scale. Value ranges are higher than SFig1
because lenient scoring raises every metric. Cells where a tool produced no
output are hatched with an en dash.
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
OUTDIR = OUT / 'SFig2'
from palette import (TOOL_ORDER, F1_CMAP, PRECISION_CMAP, RECALL_CMAP, NEUTRAL,
                     set_style, auto_text_color, save_figure)

set_style()
mpl.rcParams['hatch.color']     = '#B5B5B5'
mpl.rcParams['hatch.linewidth'] = 0.45

PROFILE = DATA / 'fusion_profile_stat_summary.tsv'
STAT    = DATA / 'stat_summary.tsv'

DEPTH_ORDER = ['1GB', '10GB', '100GB']
ID_ORDER    = ['85%', '90%', '95%']
COND_ORDER  = [(d, i) for d in DEPTH_ORDER for i in ID_ORDER]
col_keys    = [f'{d}|{i}' for d, i in COND_ORDER]


# ---------- Recompute metrics with partial calls as true recall ----------
prof = pd.read_csv(PROFILE, sep='\t')
prof = prof[prof.control == 'positive']
stat = pd.read_csv(STAT, sep='\t')
# Spiked fusion count is per depth x identity (constant across tools).
spiked = (stat[stat.control == 'positive']
          .groupby(['depth', 'Sequence_Identity'])['Spiked_Number'].first())

cat = (prof.groupby(['depth', 'Sequence_Identity', 'Algorithm', 'recall_category'])
       ['fusion_type_count'].sum().unstack(fill_value=0))
for c in ['True_Recall', 'Partial_Recall', 'False_Call']:
    if c not in cat:
        cat[c] = 0
cat = cat.reset_index()
cat['total'] = cat[['True_Recall', 'Partial_Recall', 'False_Call']].sum(axis=1)
cat['Spiked'] = cat.apply(
    lambda r: spiked.get((r['depth'], r['Sequence_Identity']), np.nan), axis=1)

tp = cat['True_Recall'] + cat['Partial_Recall']
fp = cat['False_Call']
cat['Recall'] = tp / cat['Spiked']
cat['Precision'] = np.where(tp + fp > 0, tp / (tp + fp), 0.0)
den = cat['Precision'] + cat['Recall']
cat['F1'] = np.where(den > 0, 2 * cat['Precision'] * cat['Recall'] / den, 0.0)
# A tool that produced no output at all is NA, matching SFig1.
cat.loc[cat['total'] == 0, ['Precision', 'Recall', 'F1']] = np.nan

cat['cond'] = cat.depth.astype(str) + '|' + cat.Sequence_Identity.astype(str)


def pivot(metric):
    return (cat.pivot(index='Algorithm', columns='cond', values=metric)
              .reindex(index=TOOL_ORDER, columns=col_keys))


PANELS = [
    ('Precision', pivot('Precision'), PRECISION_CMAP, 0.0, 0.9),
    ('Recall',    pivot('Recall'),    RECALL_CMAP,    0.0, 0.7),
    ('F1 score',  pivot('F1'),        F1_CMAP,        0.0, 0.7),
]


def fmt_value(v):
    """Two-decimal formatter that distinguishes a true zero from a rounded one."""
    if v == 0:
        return '0.00'
    if 0 < v < 0.005:
        return '<0.01'
    return f'{v:.2f}'


def draw_panel(fig, ax, x0, y0, mat, cell_w, cell_h, vmin, vmax, cmap,
               *, gap=0.014, font_size=8.4):
    """Paint one metric heatmap at figure-fraction origin (x0, y0)."""
    n_rows, n_cols = mat.shape
    x_pos = []
    for j, (d, _i) in enumerate(COND_ORDER):
        x_pos.append(x0 + j * cell_w + DEPTH_ORDER.index(d) * gap)

    norm = Normalize(vmin=vmin, vmax=vmax)
    for r in range(n_rows):
        y = y0 + (n_rows - 1 - r) * cell_h
        for c in range(n_cols):
            v = mat.iloc[r, c]
            x = x_pos[c]
            if pd.isna(v):
                ax.add_patch(Rectangle((x, y), cell_w, cell_h,
                             transform=fig.transFigure,
                             facecolor=NEUTRAL['na_fill'], edgecolor='white',
                             linewidth=0.9, hatch='////', zorder=2))
                ax.text(x + cell_w / 2, y + cell_h / 2, '—',
                        transform=fig.transFigure, ha='center', va='center',
                        fontsize=11, color=NEUTRAL['na_text'], zorder=3)
            else:
                ax.add_patch(Rectangle((x, y), cell_w, cell_h,
                             transform=fig.transFigure,
                             facecolor=cmap(norm(v)), edgecolor='white',
                             linewidth=0.9, zorder=2))
                ax.text(x + cell_w / 2, y + cell_h / 2, fmt_value(v),
                        transform=fig.transFigure, ha='center', va='center',
                        fontsize=font_size,
                        color=auto_text_color(v, vmin, vmax, threshold=0.50),
                        zorder=3)
    return x_pos


# ---------- Figure ----------
fig = plt.figure(figsize=(8.4, 9.6))
fig.patch.set_facecolor('white')

ax = fig.add_axes([0, 0, 1, 1])
ax.set_axis_off()
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

left   = 0.16
cell_w = 0.072
cell_h = 0.040
n_rows = len(TOOL_ORDER)
panel_h = cell_h * n_rows
panel_gap = 0.052

bottom_panel_y = 0.060
middle_panel_y = bottom_panel_y + panel_h + panel_gap
top_panel_y    = middle_panel_y + panel_h + panel_gap
panel_ys = [top_panel_y, middle_panel_y, bottom_panel_y]

x_ref = None
for (label, mat, cmap, vmin, vmax), y0 in zip(PANELS, panel_ys):
    xs = draw_panel(fig, ax, left, y0, mat, cell_w, cell_h, vmin, vmax, cmap)
    x_ref = xs

    ax.text(0.030, y0 + panel_h / 2, label, transform=fig.transFigure,
            rotation=90, ha='center', va='center',
            fontsize=12, fontweight='semibold', color=NEUTRAL['ink'])

    for r, tool in enumerate(TOOL_ORDER):
        y = y0 + (n_rows - 1 - r) * cell_h + cell_h / 2
        ax.text(left - 0.012, y, tool, transform=fig.transFigure,
                ha='right', va='center', fontsize=9.0, color=NEUTRAL['body'])

    right_edge = xs[-1] + cell_w + 0.022
    ax_cb = fig.add_axes([right_edge, y0, 0.012, panel_h])
    grad = np.linspace(0, 1, 256).reshape(-1, 1)
    ax_cb.imshow(grad, aspect='auto', cmap=cmap, origin='lower',
                 extent=[0, 1, vmin, vmax])
    ax_cb.set_xticks([])
    ax_cb.yaxis.tick_right()
    ax_cb.set_yticks(np.round(np.arange(0, vmax + 1e-9, 0.2), 1))
    ax_cb.tick_params(axis='y', length=2.5, pad=2, labelsize=7.8,
                      color=NEUTRAL['rule'])
    for spine in ax_cb.spines.values():
        spine.set_visible(False)
    ax_cb.spines['right'].set_visible(True)
    ax_cb.spines['right'].set_linewidth(0.5)
    ax_cb.spines['right'].set_color(NEUTRAL['rule'])


# ----- Shared two-tier header above the top panel -----
header_top_y = top_panel_y + panel_h
ident_y      = header_top_y + 0.008
bracket_y    = ident_y + 0.020
depth_lbl_y  = bracket_y + 0.011

for c, (_d, ident) in enumerate(COND_ORDER):
    ax.text(x_ref[c] + cell_w / 2, ident_y, ident, transform=fig.transFigure,
            ha='center', va='bottom', fontsize=8.6, color=NEUTRAL['body'])

for d in DEPTH_ORDER:
    cols = [j for j, (dd, _i) in enumerate(COND_ORDER) if dd == d]
    x_lo = x_ref[cols[0]]
    x_hi = x_ref[cols[-1]] + cell_w
    ax.plot([x_lo, x_hi], [bracket_y, bracket_y], transform=fig.transFigure,
            color=NEUTRAL['rule'], linewidth=0.9, zorder=3)
    for xe in (x_lo, x_hi):
        ax.plot([xe, xe], [bracket_y, bracket_y - 0.006],
                transform=fig.transFigure, color=NEUTRAL['rule'],
                linewidth=0.9, zorder=3)
    ax.text((x_lo + x_hi) / 2, depth_lbl_y, d, transform=fig.transFigure,
            ha='center', va='bottom', fontsize=10.0, fontweight='semibold',
            color=NEUTRAL['ink'])


# ----- Title -----
ax.text(0.045, 0.972,
        'Precision, recall and F1 with partial calls scored as true recall',
        transform=fig.transFigure, ha='left', va='top',
        fontsize=12.5, fontweight='semibold', color=NEUTRAL['ink'])


save_figure(fig, 'SFig2', OUTDIR)
plt.close(fig)
