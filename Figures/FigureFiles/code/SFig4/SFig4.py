"""Supplementary Figure 4. False, partial and true fusion calls across the simulation grid.

Three panels giving, for each of the six tools and nine simulated conditions
(three sequencing depths by three read identities), the number of false calls,
partial recalls and true recalls. Categories follow the profile classification:
false = a false-fusion type; partial = truncated, reversed, chromosomally
misaligned or wrong order; true = a correct recall. Counts are summed over fusion
subtypes. The false-call panel is a facet of two rows, the negative and positive
controls, sharing one colour scale (the negative control is a pure false-call
test). Partial and true recall come from the positive control, where alone there
is truth to recall. Each panel has its own square-root colour norm so lower
counts stay distinguishable, and cells are annotated with the exact count. Cells
where a tool produced no output are hatched with an en dash.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.colors import PowerNorm
from matplotlib.cm import ScalarMappable

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'SFig4'
from palette import (TOOL_ORDER, PRECISION_CMAP, RECALL_CMAP, F1_CMAP,
                     NEUTRAL, set_style, save_figure)

set_style()
mpl.rcParams['hatch.color']     = '#B5B5B5'
mpl.rcParams['hatch.linewidth'] = 0.45

PROFILE = DATA / 'fusion_profile_stat_summary.tsv'

DEPTH_ORDER = ['1GB', '10GB', '100GB']
ID_ORDER    = ['85%', '90%', '95%']
COND_ORDER  = [(d, i) for d in DEPTH_ORDER for i in ID_ORDER]
col_keys    = [f'{d}|{i}' for d, i in COND_ORDER]


# ---------- Data: counts per category, per tool and condition ----------
prof = pd.read_csv(PROFILE, sep='\t')
prof['cond'] = prof.depth.astype(str) + '|' + prof.Sequence_Identity.astype(str)

# Total output per (tool, condition) drives the no-output NA test: the positive
# control alone for partial/true, both controls for the pooled false calls.
pos_total = (prof[prof.control == 'positive']
             .groupby(['Algorithm', 'cond'])['fusion_type_count'].sum())
neg_total = (prof[prof.control == 'negative']
             .groupby(['Algorithm', 'cond'])['fusion_type_count'].sum())


def matrix(category, controls, na_total):
    sub = prof[(prof.recall_category == category) & prof.control.isin(controls)]
    m = (sub.groupby(['Algorithm', 'cond'])['fusion_type_count'].sum()
         .reset_index()
         .pivot(index='Algorithm', columns='cond', values='fusion_type_count')
         .reindex(index=TOOL_ORDER, columns=col_keys)).astype(float)
    m = m.where(m.notna(), 0.0)
    for tool in TOOL_ORDER:
        for ck in col_keys:
            if na_total.get((tool, ck), 0) == 0:
                m.loc[tool, ck] = np.nan
    return m


false_neg = matrix('False_Call', ['negative'], neg_total)
false_pos = matrix('False_Call', ['positive'], pos_total)
partial   = matrix('Partial_Recall', ['positive'], pos_total)
true_r    = matrix('True_Recall', ['positive'], pos_total)

# The two false-call panels share a colour scale so the negative and positive
# controls are directly comparable.
false_vmax = float(max(np.nanmax(false_neg.values), np.nanmax(false_pos.values)))

# (label, matrix, cmap, vmax, kind). The two 'false' rows form one facet panel
# (tight gap + one shared colourbar).
PANELS = [
    ('False calls\n(negative)', false_neg, PRECISION_CMAP, false_vmax, 'false'),
    ('False calls\n(positive)', false_pos, PRECISION_CMAP, false_vmax, 'false'),
    ('Partial recall', partial, RECALL_CMAP, float(np.nanmax(partial.values)), 'single'),
    ('True recall',    true_r,  F1_CMAP,     float(np.nanmax(true_r.values)),  'single'),
]


def draw_panel(fig, ax, x0, y0, mat, cell_w, cell_h, norm, cmap,
               *, gap=0.014, font_size=8.4):
    """Paint one count heatmap at figure-fraction origin (x0, y0)."""
    n_rows, n_cols = mat.shape
    x_pos = [x0 + j * cell_w + DEPTH_ORDER.index(d) * gap
             for j, (d, _i) in enumerate(COND_ORDER)]
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
                inten = norm(v)
                ax.add_patch(Rectangle((x, y), cell_w, cell_h,
                             transform=fig.transFigure,
                             facecolor=cmap(inten), edgecolor='white',
                             linewidth=0.9, zorder=2))
                ax.text(x + cell_w / 2, y + cell_h / 2, f'{int(v)}',
                        transform=fig.transFigure, ha='center', va='center',
                        fontsize=font_size, zorder=3,
                        color='white' if inten > 0.55 else NEUTRAL['ink'])
    return x_pos


# ---------- Figure ----------
fig = plt.figure(figsize=(8.4, 12.0))
fig.patch.set_facecolor('white')

ax = fig.add_axes([0, 0, 1, 1])
ax.set_axis_off()
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

left   = 0.16
cell_w = 0.072
cell_h = 0.0312
n_rows = len(TOOL_ORDER)
panel_h = cell_h * n_rows

# The two false-call rows sit in one facet (small gap); a wider gap separates the
# false, partial and true panels.
gap_in, gap_out = 0.009, 0.046
true_y      = 0.045
partial_y   = true_y + panel_h + gap_out
false_pos_y = partial_y + panel_h + gap_out
false_neg_y = false_pos_y + panel_h + gap_in
panel_ys = [false_neg_y, false_pos_y, partial_y, true_y]


def colourbar(vmax, cmap, x, y, h):
    ax_cb = fig.add_axes([x, y, 0.012, h])
    sm = ScalarMappable(norm=PowerNorm(gamma=0.5, vmin=0, vmax=vmax), cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax_cb)
    cb.outline.set_visible(False)
    cb.set_ticks(sorted({0, int(round(vmax * 0.25)), int(round(vmax * 0.5)),
                         int(round(vmax))}))
    cb.ax.tick_params(labelsize=7.8, length=2.5, pad=2, color=NEUTRAL['rule'])


x_ref = None
for (label, mat, cmap, vmax, kind), y0 in zip(PANELS, panel_ys):
    norm = PowerNorm(gamma=0.5, vmin=0, vmax=vmax)
    x_ref = draw_panel(fig, ax, left, y0, mat, cell_w, cell_h, norm, cmap)

    for r, tool in enumerate(TOOL_ORDER):
        y = y0 + (n_rows - 1 - r) * cell_h + cell_h / 2
        ax.text(left - 0.012, y, tool, transform=fig.transFigure,
                ha='right', va='center', fontsize=9.0, color=NEUTRAL['body'])

    ax.text(0.028, y0 + panel_h / 2, label, transform=fig.transFigure,
            rotation=90, ha='center', va='center',
            fontsize=9.8, fontweight='semibold', color=NEUTRAL['ink'],
            linespacing=0.95)

    # Partial and true each get a colourbar; the two false rows share one.
    if kind == 'single':
        colourbar(vmax, cmap, x_ref[-1] + cell_w + 0.022, y0, panel_h)

# One shared colourbar spanning the two-control false-call facet.
colourbar(false_vmax, PRECISION_CMAP, x_ref[-1] + cell_w + 0.022,
          false_pos_y, (false_neg_y + panel_h) - false_pos_y)


# ----- Shared two-tier header above the top panel -----
header_top_y = panel_ys[0] + panel_h
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
        'False, partial and true fusion calls across the simulation grid',
        transform=fig.transFigure, ha='left', va='top',
        fontsize=12.5, fontweight='semibold', color=NEUTRAL['ink'])


save_figure(fig, 'SFig4', OUTDIR)
plt.close(fig)
