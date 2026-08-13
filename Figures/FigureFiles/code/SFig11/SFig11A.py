"""Supplementary Figure 11A. Atypical fusion calls in Huh7.

Counts of three obvious false-chimera classes (Mitochondrial:Genomic,
Mitochondrial:Mitochondrial, Self-Misalignment) per library and replicate on a
log count axis. Columns are the four long-read libraries, rows are the three
classes and the x-axis is replicate B1 or B2. Of the long-read callers only
FusionSeeker and LongGF produce these classes, so only those two appear;
Illumina, where Arriba also reports self-misalignments, is not covered here.
Calls are kept at a minimum of two spanning reads or read pairs.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import FUSIONS, OUT
OUTDIR = OUT / 'SFig11'
from palette import TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

df = pd.read_csv(FUSIONS / 'fusions_readsupport_Huh7_discovery.tsv.gz',
                 sep='\t', low_memory=False)
reads = pd.to_numeric(df['spanning.reads'], errors='coerce').fillna(0)
pairs = pd.to_numeric(df['spanning.pairs'], errors='coerce').fillna(0)
df = df[(reads >= 2) | (pairs >= 2)]


def _canon(fid):
    """Canonical gene-pair id (A::B and B::A collapse to one); None if the id is
    missing or empty. The discovery table is not collapsed, so it carries duplicate
    rows; counting unique canonical ids rather than rows keeps duplicates from
    inflating the count."""
    if not isinstance(fid, str):
        return None
    parts = fid.replace('::', ':').split(':')
    return '::'.join(sorted(parts)) if len(parts) >= 2 and all(parts) else None


# FusionSeeker's fusion.ens.gene.id is NA for ~all its rows (an upstream bug in
# CheckKnownFusions.R sets it from a join-suffixed column); its real ENSG pair sits
# in fusion.gene.id. Coalesce so FusionSeeker is not silently dropped to zero. A no-op
# for every other caller, whose fusion.ens.gene.id is populated.
fid = df['fusion.ens.gene.id'].where(df['fusion.ens.gene.id'].notna(), df['fusion.gene.id'])
df['fusion_key'] = fid.map(_canon)
df = df.dropna(subset=['fusion_key'])

ATYPICAL_TYPES = ['Mitochondrial:Genomic',
                  'Mitochondrial:Mitochondrial',
                  'Self-Misalignment']
TOOLS = ['FusionSeeker', 'LongGF']

df = df[df.fusionType.isin(ATYPICAL_TYPES) & df.Algorithm.isin(TOOLS)]

agg = (df.groupby(['Platform', 'library_type', 'RNA_sample',
                   'Algorithm', 'fusionType'])['fusion_key']
         .nunique().reset_index(name='count'))

# Columns and rows
COLS = [
    ('PacBio', 'PCR_cDNA',   'PacBio PCR cDNA'),
    ('ONT',    'PCR_cDNA',   'ONT PCR cDNA'),
    ('ONT',    'direct_cDNA','ONT Direct cDNA'),
    ('ONT',    'direct_RNA', 'ONT Direct RNA'),
]
ROWS = ATYPICAL_TYPES
ROW_LABELS = {
    'Mitochondrial:Genomic':       'Mitochondrial:\nGenomic',
    'Mitochondrial:Mitochondrial': 'Mitochondrial:\nMitochondrial',
    'Self-Misalignment':           'Self-\nMisalignment',
}

REPS = ['B1', 'B2']


# ---------- Figure ----------
n_cols = len(COLS)
n_rows = len(ROWS)

fig = plt.figure(figsize=(11.0, 6.4))
fig.patch.set_facecolor('white')

L, R = 0.075, 0.910
B, T = 0.110, 0.770
GAP_X = 0.012
GAP_Y = 0.030
PANEL_W = (R - L - GAP_X * (n_cols - 1)) / n_cols
PANEL_H = (T - B - GAP_Y * (n_rows - 1)) / n_rows

X_POS = {'B1': 0, 'B2': 1}

axes = []
for ri, ftype in enumerate(ROWS):
    for ci, (plat, lib, _col_title) in enumerate(COLS):
        x = L + ci * (PANEL_W + GAP_X)
        y = T - (ri + 1) * PANEL_H - ri * GAP_Y
        ax = fig.add_axes([x, y, PANEL_W, PANEL_H])
        axes.append((ax, ri, ci, ftype, plat, lib))


for (ax, ri, ci, ftype, plat, lib) in axes:
    ax.set_yscale('log')
    ax.set_ylim(0.8, 5e3)
    ax.set_xlim(-0.45, 1.45)

    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(1.0)
        ax.spines[sp].set_color(NEUTRAL['ink'])
    ax.grid(True, axis='y', color=NEUTRAL['grid'], linewidth=0.3, zorder=0)

    # Bottom row: x-ticks B1/B2; others hide
    if ri == n_rows - 1:
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['B1', 'B2'], fontsize=9.4, color=NEUTRAL['body'])
    else:
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['', ''])
        ax.tick_params(axis='x', length=0)
    # Leftmost col: y-ticks; others hide
    if ci == 0:
        ax.tick_params(axis='y', labelsize=8.6, color=NEUTRAL['ink'])
    else:
        ax.tick_params(axis='y', labelleft=False, length=2)

    # Column header (top row only)
    if ri == 0:
        ax.text(0.5, 1.10, COLS[ci][2],
                transform=ax.transAxes,
                ha='center', va='bottom',
                fontsize=9.4, fontweight='semibold',
                color=NEUTRAL['ink'])

    # Row label (right column only)
    if ci == n_cols - 1:
        ax.text(1.18, 0.5, ROW_LABELS[ftype],
                transform=ax.transAxes,
                ha='left', va='center',
                fontsize=9.4, fontweight='semibold',
                color=NEUTRAL['ink'])

    # Plot points
    for rep in REPS:
        for tool in TOOLS:
            sub = agg[(agg.Platform == plat) & (agg.library_type == lib)
                      & (agg.RNA_sample == rep) & (agg.fusionType == ftype)
                      & (agg.Algorithm == tool)]
            if not len(sub):
                continue
            c = sub['count'].iloc[0]
            ax.scatter(X_POS[rep], c, s=85,
                       c=TOOL_COLORS[tool], edgecolors='black',
                       linewidths=0.7, alpha=0.92, zorder=4, clip_on=False)


# Shared axis labels
fig.text(L + (R - L) / 2, 0.045,
         'Replicate',
         ha='center', va='center',
         fontsize=10.4, fontweight='bold', color=NEUTRAL['ink'])
fig.text(0.015, B + (T - B) / 2,
         'Count  (log)',
         ha='center', va='center', rotation=90,
         fontsize=10.4, fontweight='bold', color=NEUTRAL['ink'])


# ---------- Title block + legend ----------
ax_lab = fig.add_axes([0, 0, 1, 1])
ax_lab.set_axis_off()
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)

ax_lab.text(0.035, 0.952,
            'Atypical fusion calls in Huh7',
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=12.6, fontweight='semibold', color=NEUTRAL['ink'])

# Tool color legend (only 2 tools here)
leg_y = 0.866
ax_lab.scatter(0.040, leg_y, s=90, c=TOOL_COLORS['FusionSeeker'],
               edgecolors='white', linewidths=0.7,
               transform=fig.transFigure, zorder=4, clip_on=False)
ax_lab.text(0.052, leg_y, 'FusionSeeker',
            transform=fig.transFigure,
            ha='left', va='center', fontsize=9.0, color=NEUTRAL['body'])
ax_lab.scatter(0.160, leg_y, s=90, c=TOOL_COLORS['LongGF'],
               edgecolors='white', linewidths=0.7,
               transform=fig.transFigure, zorder=4, clip_on=False)
ax_lab.text(0.172, leg_y, 'LongGF',
            transform=fig.transFigure,
            ha='left', va='center', fontsize=9.0, color=NEUTRAL['body'])


save_figure(fig, 'SFig11A', OUTDIR)
plt.close(fig)
