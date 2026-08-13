"""Figure 4A. Read length and alignment identity by library preparation.

Identical design to E2v3, but both metrics are computed over EVERY primary
alignment per library (no sampling). Each violin is a count-weighted KDE of the
pooled-replicate histogram in `huh7_readqc_hist_allreads.tsv`; medians are exact
weighted medians. ONT blue / PacBio red. Illumina absent.
"""
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig4'
from palette import (NEUTRAL, set_style, save_figure, PLATFORM_COLORS, PREP_ORDER,
                     prep_key)

set_style()

HIST_PATH = DATA / 'huh7_readqc_hist_allreads.tsv'

h = pd.read_csv(HIST_PATH, sep='\t')
h['prep'] = [prep_key(p, l) for p, l in zip(h['Platform'], h['library_type'])]
g = h.groupby(['prep', 'metric', 'bin'], as_index=False)['count'].sum()  # pool reps

PREPS = [p for p in PREP_ORDER if p in set(g['prep'])]


def platform_of(prep):
    return prep.split()[0]


def wmedian(values, weights):
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    o = np.argsort(values)
    v, cw = values[o], np.cumsum(weights[o])
    return float(v[np.searchsorted(cw, cw[-1] / 2.0)])


def series(prep, metric):
    s = g[(g.prep == prep) & (g.metric == metric)]
    return s['bin'].to_numpy(float), s['count'].to_numpy(float)


Y = {p: len(PREPS) - 1 - i for i, p in enumerate(PREPS)}
HALFW = 0.38
GRID = 400


def draw_violin(ax, xgrid, xplot, kde_values, kde_weights, ybase, colour,
                median_x, median_kde_x, label, label_fmt):
    kde = gaussian_kde(kde_values, weights=kde_weights)
    dmax = float(kde(xgrid).max())
    dens = kde(xgrid) / dmax * HALFW
    ax.fill_between(xplot, ybase - dens, ybase + dens, facecolor=colour,
                    edgecolor=colour, linewidth=0.8, alpha=0.85, zorder=2)
    hmed = float(kde(median_kde_x)) / dmax * HALFW
    ax.plot([median_x, median_x], [ybase - hmed, ybase + hmed], color=NEUTRAL['ink'],
            linewidth=1.5, zorder=4,
            path_effects=[pe.withStroke(linewidth=3.0, foreground='white')])
    ax.text(median_x, ybase + HALFW + 0.06, label_fmt.format(label), ha='center',
            va='bottom', fontsize=8.3, fontweight='semibold', color=NEUTRAL['ink'],
            zorder=5)


fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

LEN_XLIM = (90, 20000)
lg_grid = np.linspace(np.log10(LEN_XLIM[0] * 0.95), np.log10(LEN_XLIM[1] * 1.05), GRID)
len_xplot = 10 ** lg_grid
for prep in PREPS:
    vals, cnts = series(prep, 'length')
    med = wmedian(vals, cnts)
    draw_violin(axL, lg_grid, len_xplot, np.log10(vals), cnts, Y[prep],
                PLATFORM_COLORS[platform_of(prep)],
                median_x=med, median_kde_x=np.log10(med),
                label=med, label_fmt='{:.0f} bp')

axL.set_xscale('log')
axL.set_xlim(*LEN_XLIM)
axL.set_xticks([100, 300, 1000, 3000, 10000])
axL.set_xticklabels(['100', '300', '1,000', '3,000', '10,000'])
axL.set_xlabel('Read length (bp), log scale')

ID_XLIM = (83, 101.3)
id_grid = np.linspace(83.4, 100.35, GRID)
for prep in PREPS:
    vals, cnts = series(prep, 'identity')
    vals = vals / 10.0
    med = wmedian(vals, cnts)
    draw_violin(axR, id_grid, id_grid, vals, cnts, Y[prep],
                PLATFORM_COLORS[platform_of(prep)],
                median_x=med, median_kde_x=med,
                label=med, label_fmt='{:.1f}%')

axR.set_xlim(*ID_XLIM)
axR.set_xticks([85, 90, 95, 100])
axR.set_xlabel('Per-read alignment identity (%)')

for ax in (axL, axR):
    ax.set_ylim(-0.6, len(PREPS) - 0.4)
    ax.grid(axis='x', color=NEUTRAL['grid'], linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis='y', length=0)

axL.set_yticks([Y[p] for p in PREPS])
axL.set_yticklabels(PREPS)
for tick in axL.get_yticklabels():
    tick.set_color(NEUTRAL['ink'])
    tick.set_fontweight('bold')
    tick.set_fontsize(9.5)
plt.setp(axR.get_yticklabels(), visible=False)

fig.suptitle('Read length and identity by library preparation',
             x=0.012, y=0.985, ha='left', fontsize=12.8, fontweight='semibold',
             color=NEUTRAL['ink'])

_plat = ['ONT', 'PacBio']
_ph = [Line2D([0], [0], marker='o', linestyle='none', markersize=8,
              markerfacecolor=PLATFORM_COLORS[p], markeredgecolor='white',
              markeredgewidth=0.7, label=p) for p in _plat]
_leg = fig.legend(handles=_ph, loc='upper right', bbox_to_anchor=(0.985, 0.998),
                  frameon=False, ncol=2, handletextpad=0.35, columnspacing=1.1,
                  fontsize=9.0)
for _t, _p in zip(_leg.get_texts(), _plat):
    _t.set_color(PLATFORM_COLORS[_p])
    _t.set_fontweight('bold')

fig.subplots_adjust(left=0.135, right=0.985, top=0.9, bottom=0.115, wspace=0.08)
save_figure(fig, 'Fig4A', OUTDIR)
print('E2v4 rendered')
