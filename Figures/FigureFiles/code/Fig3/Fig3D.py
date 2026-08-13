#!/usr/bin/env python3
"""Figure 3D. Distinct-fusion burden by library preparation.

Strip plot of distinct fusions per library (a fusion = a distinct gene pair, kept
at spanning.reads >= 2 or spanning.pairs >= 2) for each caller, grouped along x by
library preparation in escalating order of amplification: ONT direct RNA, ONT
direct cDNA, ONT PCR-cDNA, PacBio PCR-cDNA and finally Illumina PCR-cDNA. One point
per caller x replicate, coloured by caller with a deterministic dodge. A black
crossbar marks the per-prep median and a faint line traces the trend. The point:
the number of distinct fusions balloons with amplification and short reads, running
from tens to hundreds of thousands.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from matplotlib.ticker import LogLocator, FuncFormatter, NullFormatter

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig3'
from palette import TOOL_COLORS, SHORT_READ_COLORS, NEUTRAL, set_style, save_figure

set_style()
ALLC = {**TOOL_COLORS, **SHORT_READ_COLORS}

df = pd.read_csv(DATA / 'call_burden.tsv', sep='\t')
df['prep'] = list(zip(df['Platform'], df['library_type']))

# Preparations in escalating order of amplification / read shortening.
preps = [
    ('ONT', 'direct_RNA'),
    ('ONT', 'direct_cDNA'),
    ('ONT', 'PCR_cDNA'),
    ('PacBio', 'PCR_cDNA'),
    ('Illumina', 'PCR_cDNA'),
]
prep_labels = [
    'ONT\ndirect RNA',
    'ONT\ndirect cDNA',
    'ONT\nPCR-cDNA',
    'PacBio\nPCR-cDNA',
    'Illumina\nPCR-cDNA',
]
prep_ix = {p: i for i, p in enumerate(preps)}
assert df['prep'].map(prep_ix).notna().all(), 'unmapped preparation in data'

caller_order = ['CTAT-LR-Fusion', 'JAFFAL', 'Genion', 'FusionSeeker', 'GFSeeker', 'LongGF']
short_read = ['JAFFA-direct', 'Arriba', 'STAR-Fusion']   # Illumina-only comparators
legend_order = [t for t in caller_order + short_read if t in set(df['Algorithm'])]

# Deterministic horizontal dodge: callers spread within each prep, the two
# replicates sit either side of the caller centre. The short-read tools appear
# only in the Illumina column, so they spread within it.
dodge = dict(zip(caller_order, np.linspace(-0.26, 0.26, len(caller_order))))
for t, dx in zip(short_read, np.linspace(-0.15, 0.15, len(short_read))):
    dodge[t] = dx
rep_dx = {'B1': -0.03, 'B2': 0.03}

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.set_axisbelow(True)

# Faint platform boundaries (ONT | PacBio | Illumina).
for xb in (2.5, 3.5):
    ax.axvline(xb, color=NEUTRAL['grid'], lw=0.8, zorder=0)

# One point per caller x replicate.
for _, r in df.iterrows():
    x = prep_ix[r['prep']] + dodge[r['Algorithm']] + rep_dx[r['RNA_sample']]
    ax.scatter(x, r['n_fusions'], s=50, color=ALLC[r['Algorithm']],
               edgecolor='white', linewidth=0.5, alpha=0.95, zorder=7)

# Per-prep median crossbar and the trend that connects them.
medians = [float(np.median(df[df['prep'] == p]['n_fusions'])) for p in preps]
ax.plot(range(4), medians[:4], color=NEUTRAL['muted'], lw=1.2, alpha=0.85, zorder=3)
ax.plot([3, 4], medians[3:], color=NEUTRAL['muted'], lw=1.2, ls=(0, (4, 3)),
        alpha=0.7, zorder=3)
for i, med in enumerate(medians):
    ax.plot([i - 0.34, i + 0.34], [med, med], color=NEUTRAL['ink'], lw=2.3,
            solid_capstyle='round', zorder=6,
            path_effects=[pe.withStroke(linewidth=4.4, foreground='white')])

# ----- axes -----
ax.set_yscale('log')
ax.set_ylim(4, 1.5e6)
ax.set_xlim(-0.55, 4.55)
ax.yaxis.set_major_locator(LogLocator(base=10))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{int(v):,}'))
ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(2, 3, 4, 5, 6, 7, 8, 9), numticks=12))
ax.yaxis.set_minor_formatter(NullFormatter())
ax.tick_params(axis='y', which='minor', length=2, color=NEUTRAL['rule'])
ax.grid(axis='y', which='major', color=NEUTRAL['grid'], lw=0.7, zorder=0)

ax.set_xticks(range(5))
ax.set_xticklabels(prep_labels, linespacing=1.25)
ax.tick_params(axis='x', length=0, pad=6)
ax.set_ylabel('Distinct fusions per library (log)')

ax.set_title('Distinct fusions per library balloon with amplification',
             loc='left', fontsize=12.8, fontweight='semibold',
             color=NEUTRAL['ink'], pad=10)

# ----- legend: tool names bold and in the tool colour -----
handles = [Line2D([0], [0], marker='o', linestyle='none', markersize=7,
                  markerfacecolor=ALLC[t], markeredgecolor='white',
                  markeredgewidth=0.5, label=t) for t in legend_order]
leg = ax.legend(handles=handles, loc='upper left', frameon=False,
                handletextpad=0.4, labelspacing=0.52, borderaxespad=0.5,
                fontsize=8.5)
for txt, t in zip(leg.get_texts(), legend_order):
    txt.set_color(ALLC[t])
    txt.set_fontweight('bold')

print('per-prep medians:', {lab.replace(chr(10), ' '): round(m) for lab, m in zip(prep_labels, medians)})
save_figure(fig, 'Fig3D', OUTDIR)
print('done')
