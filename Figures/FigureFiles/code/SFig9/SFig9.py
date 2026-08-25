"""Supplementary Figure 9. Huh7 fusion-discovery summary per library and metric.

Rows are discovery metrics (total, known-DB matches, in Huh7 literature,
RT-PCR/Sanger validated, HCC-only literature); columns are the five
library-by-platform combinations; each cell is the count of unique fusions on a
log colour scale. The right-hand bar gives the count in the union across all
libraries (deduplicated).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import OUT
OUTDIR = OUT / 'SFig9'
from palette import F1_CMAP, NEUTRAL, set_style, save_figure
from _data import load, LIB_ORDER, KNOWN_TYPES, _fusion_key as _canon

set_style()

df = load()


def _fids(sub):
    """Unique canonical fusion ids in a slice. Gene pairs are sorted so tools
    that disagree on gene1/gene2 order collapse to one id."""
    ids = set()
    for v in sub['fusion.ens.gene.id'].dropna():
        c = _canon(v)
        if c is not None:
            ids.add(c)
    return ids


# Compute per-library values + the union across libraries
per_lib = {}
union_total, union_known, union_rtpcr = set(), set(), set()
union_huh7, union_hcc = set(), set()

for plat, lib, label in LIB_ORDER:
    sub = df[(df.Platform == plat) & (df.library_type == lib)]
    sub_k = sub[sub.Discovery.isin(KNOWN_TYPES)]
    sub_rt = sub_k[sub_k.Validation.fillna('').str.contains('RT-PCR')]
    sub_huh = sub_k[sub_k.Sample == 'HUH7']
    sub_hcc = sub_k[sub_k.Sample.fillna('').isin(
        ['HCC', 'HEPG2', 'Fibrolamellar-HCC'])]
    total_ids = _fids(sub)
    known_ids = _fids(sub_k)
    rt_ids    = _fids(sub_rt)
    huh_ids   = _fids(sub_huh)
    hcc_ids   = _fids(sub_hcc) - huh_ids   # HCC-only

    per_lib[label] = {
        'total':    len(total_ids),
        'known':    len(known_ids),
        'rtpcr':    len(rt_ids),
        'huh7_lit': len(huh_ids),
        'hcc_only': len(hcc_ids),
    }
    union_total |= total_ids
    union_known |= known_ids
    union_rtpcr |= rt_ids
    union_huh7  |= huh_ids
    union_hcc   |= hcc_ids

# Union across libraries (deduplicated)
# HCC-only at union level = HCC literature matches not in Huh7 literature
union_hcc_only = union_hcc - union_huh7

# Metrics (rows of the grid)
METRICS = [
    ('total',    'Total fusions called',     union_total),
    ('known',    'Known DB matches',         union_known),
    ('huh7_lit', 'In Huh7 literature',       union_huh7),
    ('rtpcr',    'RT-PCR / Sanger validated', union_rtpcr),
    ('hcc_only', 'HCC-only literature\n(new in Huh7)', union_hcc_only),
]

# Build matrix
LIB_LABELS = [l for _, _, l in LIB_ORDER]
M = np.zeros((len(METRICS), len(LIB_LABELS)))
for i, (key, _, _) in enumerate(METRICS):
    for j, lab in enumerate(LIB_LABELS):
        M[i, j] = per_lib[lab][key]

union_totals = [len(METRICS[i][2]) for i in range(len(METRICS))]

# Shared log color scale across the heatmap
vmax = max(M.max(), 1)
vmin = 1
norm = LogNorm(vmin=vmin, vmax=vmax)


# ---------- Figure ----------
# Heatmap cell side ~ 0.85 inch (square)
n_rows, n_cols = M.shape
CELL_INCH = 0.85

FIG_W = 13.6
FIG_H = 7.0
fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor('white')

# Layout:
#   left margin (row labels)   |   heatmap (HM)   |   small gap   |   summary bar (BAR)
ROW_LBL_W = 0.20
GAP_HM_BAR = 0.025
BAR_W_INCH = 1.6
HM_W_IN = n_cols * CELL_INCH
HM_H_IN = n_rows * CELL_INCH

# Figure-coord positions
hm_left = ROW_LBL_W
hm_w = HM_W_IN / FIG_W
hm_h = HM_H_IN / FIG_H
bar_left = hm_left + hm_w + GAP_HM_BAR
bar_w   = BAR_W_INCH / FIG_W

# Vertical: leave room below for the column labels and the colourbar
BOT = 0.220
TOP_AFTER = BOT + hm_h    # heatmap top

ax_hm = fig.add_axes([hm_left, BOT, hm_w, hm_h])

# Heatmap
im = ax_hm.imshow(np.ma.masked_where(M == 0, M),
                  cmap=F1_CMAP, norm=norm,
                  aspect='equal', interpolation='nearest', origin='upper')

# Column labels (rotated below)
ax_hm.set_xticks(range(n_cols))
ax_hm.set_xticklabels([l.replace('\n', ' ') for l in LIB_LABELS],
                      rotation=24, ha='right',
                      fontsize=10.0, color=NEUTRAL['body'])
# Row labels (left)
ax_hm.set_yticks(range(n_rows))
ax_hm.set_yticklabels([m[1] for m in METRICS],
                      fontsize=10.2, color=NEUTRAL['ink'])

# Annotations (counts)
for i in range(n_rows):
    for j in range(n_cols):
        v = M[i, j]
        if v == 0:
            txt, color, weight = '0', '#A8A8A8', 'regular'
        else:
            if v < 1000:   txt = f'{int(v)}'
            elif v < 1e4:  txt = f'{v/1000:.1f}k'
            elif v < 1e6:  txt = f'{int(v/1000)}k'
            else:          txt = f'{v/1e6:.2f}M'
            color = 'white' if norm(v) > 0.55 else NEUTRAL['ink']
            weight = 'semibold'
        ax_hm.text(j, i, txt, ha='center', va='center',
                   fontsize=10.2, fontweight=weight, color=color)

ax_hm.set_xticks(np.arange(n_cols + 1) - 0.5, minor=True)
ax_hm.set_yticks(np.arange(n_rows + 1) - 0.5, minor=True)
ax_hm.grid(which='minor', color='white', linewidth=1.0)
ax_hm.tick_params(which='minor', length=0)
ax_hm.tick_params(which='major', length=0)
for sp in ax_hm.spines.values():
    sp.set_visible(False)

# --- Right summary bar (across-library deduplicated count per metric) ---
ax_bar = fig.add_axes([bar_left, BOT, bar_w, hm_h])
ys = np.arange(n_rows)
xmax = max(union_totals)
ax_bar.set_xlim(0, xmax * 1.18)
ax_bar.set_ylim(-0.5, n_rows - 0.5)
ax_bar.invert_yaxis()

for sp in ('top', 'right', 'left'):
    ax_bar.spines[sp].set_visible(False)
ax_bar.spines['bottom'].set_linewidth(1.0)
ax_bar.spines['bottom'].set_color(NEUTRAL['ink'])
ax_bar.tick_params(axis='y', length=0, labelleft=False)
ax_bar.tick_params(axis='x', labelsize=8.6, color=NEUTRAL['ink'])
ax_bar.set_xscale('log')
ax_bar.set_xlim(0.7, xmax * 3)
ax_bar.grid(True, axis='x', color=NEUTRAL['grid'], linewidth=0.4, zorder=0)

for i, v in enumerate(union_totals):
    color = F1_CMAP(norm(max(v, 1)))
    ax_bar.barh(i, v, height=0.62, color=color,
                edgecolor=NEUTRAL['ink'], linewidth=0.6, zorder=3)
    ax_bar.text(v * 1.15, i, f'{int(v):,}',
                ha='left', va='center',
                fontsize=9.2, fontweight='bold', color=NEUTRAL['ink'])

# Header for the right-side summary bar, placed just above ax_bar.
fig.text(bar_left + bar_w / 2, BOT + hm_h + 0.018,
         'Union across\nall libraries',
         ha='center', va='bottom',
         fontsize=9.6, fontweight='semibold', color=NEUTRAL['ink'])
ax_bar.set_xlabel('Unique fusions (log)',
                  fontsize=9.0, fontweight='bold',
                  color=NEUTRAL['ink'], labelpad=4)


# ---------- Colourbar (under heatmap) ----------
ax_cb = fig.add_axes([hm_left, BOT - 0.16, hm_w, 0.018])
cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8.6, color=NEUTRAL['ink'], length=3)
cb.set_label('Per-cell fusion count (log)',
             fontsize=9.2, fontweight='bold',
             color=NEUTRAL['ink'], labelpad=4)


# ---------- Title ----------
fig.text(0.035, 0.945,
         'Huh7 fusion discovery summary per library and metric',
         ha='left', va='top',
         fontsize=12.6, fontweight='semibold', color=NEUTRAL['ink'])


save_figure(fig, 'SFig9', OUTDIR)
plt.close(fig)
