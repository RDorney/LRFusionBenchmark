"""Figure 4B. Transcript coverage completeness by library preparation.

Meta-transcript 5' to 3' coverage from the Huh7 transcriptome alignments, faceted
by transcript-length bin. Within each bin the mean per-transcript normalised
coverage (each transcript with at least 10 reads contributes its own coverage
profile with equal weight, so high-expression transcripts do not dominate) is
drawn across transcript percentiles, 0% = 5' end, 100% = 3' end. Coverage falls off with
transcript length and differs by preparation, showing which libraries deliver reads
that span whole transcripts, and hence fusion breakpoints. Colour follows the
platform identity used across the paper (ONT a blue family, PacBio red), matching
the other library-preparation panels.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig4'
from palette import NEUTRAL, set_style, save_figure

set_style()

d = pd.read_csv(DATA / 'huh7_txcoverage_pertx.tsv', sep='\t')
d['prep'] = d['Platform'] + ' ' + d['library_type'].str.replace('_', ' ')
LABELS = {'ONT PCR cDNA': 'ONT PCR-cDNA', 'PacBio PCR cDNA': 'PacBio PCR-cDNA'}
d['prep'] = d['prep'].replace(LABELS)

# Pool the two replicates: transcript-weighted mean of the per-transcript
# averaged coverage.
d['w'] = d['mean_cov'] * d['n_transcripts']
pool = (d.groupby(['prep', 'length_bin', 'percentile'])
          .agg(w=('w', 'sum'), n=('n_transcripts', 'sum')).reset_index())
pool['norm'] = pool['w'] / pool['n']

# ONT chemistries as a blue family (gentlest = darkest), PacBio red.
PREP_ORDER = ['ONT direct RNA', 'ONT direct cDNA', 'ONT PCR-cDNA', 'PacBio PCR-cDNA']
PREP_COLOR = {
    'ONT direct RNA':  '#0D3B54',
    'ONT direct cDNA': '#2C7FB8',
    'ONT PCR-cDNA':    '#7FBEE0',
    'PacBio PCR-cDNA': '#C73E5C',
}
BINS = ['0-1kb', '1-2.5kb', '2.5kb+']
BIN_TITLES = {'0-1kb': '0-1 kb', '1-2.5kb': '1-2.5 kb', '2.5kb+': '2.5 kb+'}

fig, axes = plt.subplots(1, 3, figsize=(11, 4.3), sharey=True)
for ax, b in zip(axes, BINS):
    for prep in PREP_ORDER:
        g = pool[(pool.prep == prep) & (pool.length_bin == b)].sort_values('percentile')
        if len(g):
            ax.plot(g['percentile'], g['norm'], color=PREP_COLOR[prep], lw=2.1,
                    solid_capstyle='round', zorder=3)
    ax.set_title(BIN_TITLES[b], fontsize=10.6, fontweight='semibold',
                 color=NEUTRAL['ink'], pad=6)
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.02, 1.04)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.grid(True, color=NEUTRAL['grid'], lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel("Position along transcript (%, 5' to 3')")

axes[0].set_ylabel('Normalised coverage')

# Centred group label for the facets.
fig.text(0.525, 0.895, 'Transcript length', ha='center', va='bottom',
         fontsize=9.4, color=NEUTRAL['muted'])

# Legend: platform-coloured, bold prep names.
handles = [Line2D([0], [0], color=PREP_COLOR[p], lw=2.6, label=p) for p in PREP_ORDER]
leg = fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.992, 0.998),
                 frameon=False, ncol=2, fontsize=8.6, handlelength=1.7,
                 columnspacing=1.3, handletextpad=0.5)
for t, p in zip(leg.get_texts(), PREP_ORDER):
    t.set_color(PREP_COLOR[p])
    t.set_fontweight('bold')

fig.suptitle('Transcript coverage completeness by library preparation',
             x=0.012, y=0.985, ha='left', fontsize=12.8, fontweight='semibold',
             color=NEUTRAL['ink'])
fig.subplots_adjust(left=0.065, right=0.985, top=0.80, bottom=0.14, wspace=0.09)
save_figure(fig, 'Fig4B', OUTDIR)
print('done')
