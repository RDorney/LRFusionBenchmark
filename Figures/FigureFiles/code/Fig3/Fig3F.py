"""Figure 3F. Fusion-partner promiscuity by caller and library preparation.

Fusion-promiscuity dumbbell: the typical gene stays paired, a few hub genes run away.

One dumbbell per caller x library-prep (replicates averaged): a small open dot at
the MEAN partners per gene (the typical gene) linked to a filled dot at the MAX
(the worst hub gene). Preps run top-to-bottom from gentlest to harshest. In every
long-read prep the mean sits at ~1 partner while the max runs to hundreds — a few
rogue hubs on an otherwise healthy library. Only Illumina short-read pushes the
whole library up: even its typical gene reaches ~100 partners.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig3'
from palette import TOOL_COLORS, SHORT_READ_COLORS, NEUTRAL, set_style, save_figure

set_style()
ALLC = {**TOOL_COLORS, **SHORT_READ_COLORS}

d = pd.read_csv(DATA / 'promiscuity_lib_metrics.tsv', sep='\t')
agg = (d.groupby(['prep', 'Algorithm'])
         .agg(mean_p=('mean_partners', 'mean'), max_p=('max_partners', 'mean'))
         .reset_index())

PREPS = ['ONT direct RNA', 'ONT direct cDNA', 'ONT PCR cDNA',
         'PacBio PCR cDNA', 'Illumina PCR cDNA']
PREP_LAB = {'ONT direct RNA': 'ONT\ndirect RNA',
            'ONT direct cDNA': 'ONT\ndirect cDNA',
            'ONT PCR cDNA': 'ONT\nPCR-cDNA',
            'PacBio PCR cDNA': 'PacBio\nPCR-cDNA',
            'Illumina PCR cDNA': 'Illumina\nPCR-cDNA'}
TOOLS = ['STAR-Fusion', 'LongGF', 'CTAT-LR-Fusion', 'Genion', 'Arriba', 'JAFFAL',
         'GFSeeker', 'FusionSeeker', 'JAFFA-direct']   # cleanest -> noisiest (top -> bottom)

present = {p: [t for t in TOOLS if ((agg.prep == p) & (agg.Algorithm == t)).any()]
           for p in PREPS}
heights = [max(len(present[p]), 2.4) for p in PREPS]

XMIN, XMAX = 0.8, 22000
fig = plt.figure(figsize=(10.6, 9.0))
gs = fig.add_gridspec(len(PREPS), 1, left=0.185, right=0.95, top=0.845,
                      bottom=0.085, hspace=0.22, height_ratios=heights)

axes = []
for gi, prep in enumerate(PREPS):
    ax = fig.add_subplot(gs[gi, 0], sharex=axes[0] if axes else None)
    axes.append(ax)
    tools = present[prep]
    n = len(tools)
    ys = list(range(n))[::-1]                 # first tool at the top
    ax.set_ylim(-0.75, n - 0.25)
    ax.set_xscale('log')
    ax.set_xlim(XMIN, XMAX)
    ax.axvline(1, color=NEUTRAL['rule'], lw=1.0, linestyle=(0, (4, 3)), zorder=1)

    for t, y in zip(tools, ys):
        r = agg[(agg.prep == prep) & (agg.Algorithm == t)].iloc[0]
        c = ALLC[t]
        mp, xp = float(r.mean_p), float(r.max_p)
        ax.plot([mp, xp], [y, y], color=c, lw=2.6, alpha=0.5,
                solid_capstyle='round', zorder=2)
        ax.scatter([mp], [y], s=36, facecolor='white', edgecolor=c,
                   linewidths=1.5, zorder=3)
        ax.scatter([xp], [y], s=128, facecolor=c, edgecolor='white',
                   linewidths=0.8, zorder=4)
        ax.annotate(f"{xp:,.0f}", (xp, y), textcoords='offset points',
                    xytext=(9, 0), va='center', ha='left', fontsize=8.6,
                    fontweight='semibold', color=c, zorder=5)

    ax.set_yticks(ys)
    for tl, t in zip(ax.set_yticklabels(tools, fontsize=9.3, fontweight='bold'), tools):
        tl.set_color(ALLC[t])
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color(NEUTRAL['ink'])
        ax.spines[sp].set_linewidth(1.0)
    ax.tick_params(axis='y', length=0)
    ax.grid(True, axis='x', color=NEUTRAL['grid'], lw=0.4, zorder=0)
    if gi < len(PREPS) - 1:
        ax.tick_params(axis='x', labelbottom=False)

# shared x ticks on the bottom axis
axes[-1].xaxis.set_major_locator(plt.FixedLocator([1, 10, 100, 1000, 10000]))
axes[-1].xaxis.set_major_formatter(plt.FixedFormatter(['1', '10', '100', '1,000', '10,000']))
axes[-1].tick_params(axis='x', labelsize=9.4, color=NEUTRAL['ink'])

# ---- vertical prep labels in the left margin ----
lab = fig.add_axes([0, 0, 1, 1]); lab.set_axis_off()
lab.set_xlim(0, 1); lab.set_ylim(0, 1)
for ax, prep in zip(axes, PREPS):
    p = ax.get_position()
    lab.text(0.072, (p.y0 + p.y1) / 2, PREP_LAB[prep], ha='center', va='center',
             rotation=90, fontsize=9.0, fontweight='semibold', color=NEUTRAL['ink'],
             linespacing=1.5)
# ---- title + key ----
lab.text(0.035, 0.965, 'The typical gene keeps one partner; a few hubs run away',
         ha='left', va='top', fontsize=13, fontweight='semibold', color=NEUTRAL['ink'])
kx = 0.275
lab.scatter(kx, 0.905, s=40, facecolor='white', edgecolor=NEUTRAL['ink'],
            linewidths=1.5, zorder=5)
lab.text(kx + 0.013, 0.905, 'mean (typical gene)', va='center', fontsize=9,
         color=NEUTRAL['body'])
lab.scatter(kx + 0.235, 0.905, s=120, facecolor=NEUTRAL['muted'],
            edgecolor='white', linewidths=0.8, zorder=5)
lab.text(kx + 0.248, 0.905, 'max (worst hub gene)', va='center', fontsize=9,
         color=NEUTRAL['body'])

fig.text(0.61, 0.028, 'Partners per gene  (log scale)', ha='center', va='center',
         fontsize=10.6, fontweight='bold', color=NEUTRAL['ink'])

save_figure(fig, 'Fig3F', OUTDIR)
print('done')
