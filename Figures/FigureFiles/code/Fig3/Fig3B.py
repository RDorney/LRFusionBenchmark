"""Figure 3B. Novel-call heatmap by tool and fusion subtype, per cell line.

Each row is a tool and each column a fusion subtype; the cell value is the number
of unique novel fusions, deduplicated across the three library preps (full-depth
SGNex data), on a log colour scale. The two panels (K562, MCF7) share one colour
scale and each cell prints its count.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig3'
from palette import TOOL_ORDER, F1_CMAP, NEUTRAL, set_style, save_figure

set_style()

df = pd.read_csv(DATA / 'fusionlevel_SGNex_summary.tsv.gz', sep='\t')
df = df[df.Sequencing_Depth == 'Total'].copy()

CELL_LINES = ['K562', 'MCF7']
SUBTYPES = [
    'SAGe',
    'inter-chromosomal',
    'intra-chromosomal',
    'Mitochondrial:Genomic',
    'Mitochondrial:Mitochondrial',
    'read-through',
    'Self-Misalignment',
    'Sense-Antisense',
    'tri-fusion',
]

# Unique novel fusions across libraries: distinct fusionGeneID per (cell, tool,
# subtype) among the novel calls at full (Total) depth.
agg = (df[df.novelnumber > 0]
         .groupby(['Cell_Lines', 'Algorithm', 'fusionType'])['fusionGeneID']
         .nunique().reset_index(name='novelnumber'))

# Build a wide matrix per cell line
matrices = {}
for cell in CELL_LINES:
    m = np.zeros((len(TOOL_ORDER), len(SUBTYPES)))
    sub = agg[agg.Cell_Lines == cell]
    for _, r in sub.iterrows():
        try:
            i = TOOL_ORDER.index(r['Algorithm'])
            j = SUBTYPES.index(r['fusionType'])
        except ValueError:
            continue
        m[i, j] = r['novelnumber']
    matrices[cell] = m

# Shared log colour scale using F1_CMAP (cream to teal to navy).
all_vals = np.concatenate([matrices[c].ravel() for c in CELL_LINES])
vmax = max(int(all_vals.max()), 1)
vmin = 1.0
norm = LogNorm(vmin=vmin, vmax=vmax)
cmap = F1_CMAP


# ---------- Figure ----------
fig = plt.figure(figsize=(13.4, 5.4))
fig.patch.set_facecolor('white')

# Layout: 2 panels + a thin vertical colourbar on the far right.
# Bigger bottom margin so rotated subtype labels don't clash with anything.
L, R = 0.110, 0.900
B, T = 0.255, 0.735
GAP  = 0.060
PANEL_W = (R - L - GAP) / 2
PANEL_H = T - B


def draw_heatmap(ax, mat, cell, show_ylabel):
    nr, nc = mat.shape
    # Mask zeros (so they render as white/empty, not the lowest log colour)
    mat_masked = np.ma.masked_where(mat == 0, mat)

    im = ax.imshow(mat_masked, aspect='auto', cmap=cmap, norm=norm,
                   origin='upper', interpolation='nearest')

    ax.set_xticks(range(nc))
    ax.set_xticklabels(SUBTYPES, rotation=38, ha='right',
                       fontsize=9.0, color=NEUTRAL['body'])
    if show_ylabel:
        ax.set_yticks(range(nr))
        ax.set_yticklabels(TOOL_ORDER, fontsize=9.5, color=NEUTRAL['ink'])
    else:
        ax.set_yticks(range(nr))
        ax.set_yticklabels([''] * nr)

    # Cell value annotations
    for i in range(nr):
        for j in range(nc):
            v = mat[i, j]
            if v == 0:
                txt = '0'
                color = '#A8A8A8'
                weight = 'regular'
            else:
                if v < 1000:
                    txt = f'{int(v)}'
                elif v < 10000:
                    txt = f'{v/1000:.1f}k'
                else:
                    txt = f'{int(v/1000)}k'
                bg_intensity = norm(v)
                # F1_CMAP gets dark (navy) above ~0.55 → switch text to white
                color = 'white' if bg_intensity > 0.55 else NEUTRAL['ink']
                weight = 'semibold'
            ax.text(j, i, txt,
                    ha='center', va='center',
                    fontsize=8.6, fontweight=weight, color=color)

    # Light grid
    ax.set_xticks(np.arange(nc + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(nr + 1) - 0.5, minor=True)
    ax.grid(which='minor', color='white', linewidth=1.0)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=0)

    for sp in ax.spines.values():
        sp.set_visible(False)

    ax.set_title(cell, fontsize=12.4, fontweight='semibold',
                 color=NEUTRAL['ink'], loc='left', pad=10)

    return im


axes = []
for i, cell in enumerate(CELL_LINES):
    ax = fig.add_axes([L + i * (PANEL_W + GAP), B, PANEL_W, PANEL_H])
    im = draw_heatmap(ax, matrices[cell], cell, show_ylabel=(i == 0))
    axes.append(ax)


# ---------- Vertical colourbar on the right ----------
cbar_ax = fig.add_axes([R + 0.012, B, 0.016, PANEL_H])
cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8.6, color=NEUTRAL['ink'], length=3)
cb.set_label('Unique novel fusion calls  (count, log)',
             fontsize=9.6, fontweight='bold', color=NEUTRAL['ink'],
             labelpad=8, rotation=270, va='bottom')


# ---------- Title block ----------
ax_lab = fig.add_axes([0, 0, 1, 1])
ax_lab.set_axis_off()
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)

ax_lab.text(0.035, 0.960,
            'Novel fusion calls by tool and fusion subtype',
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])


save_figure(fig, 'Fig3B', OUTDIR)
plt.close(fig)
