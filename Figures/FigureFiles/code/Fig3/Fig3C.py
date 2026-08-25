"""Figure 3C. Huh7 fusion-call distributions across libraries.

Three stacked panels: total calls, known-database recovery, and putative novel.
In each panel the x-axis is the five Huh7 libraries and the box shows the spread
of call counts across tools and replicates, with individual tool values overlaid
as dots. Total and novel use a log scale, known recovery a linear scale. Calls
are filtered to a minimum read support of two.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import OUT
OUTDIR = OUT / 'Fig3'
from palette import TOOL_COLORS, SHORT_READ_COLORS, NEUTRAL, set_style, save_figure
from _data import load, LIB_ORDER, KNOWN_TYPES

set_style()

df = load()  # minimum read support of two

KEYS = ['Platform', 'library_type', 'Algorithm', 'RNA_sample']

# The Illumina short-read callers (JAFFA-direct, Arriba, STAR-Fusion) sit outside the
# long-read palette and carry their own tints. Without them they fell back to grey,
# which reads as LongGF (also grey) in the Illumina column.
ALL_TOOL_COLORS = {**TOOL_COLORS, **SHORT_READ_COLORS}


def _oriented(fid):
    """Oriented gene-pair key, 5' partner first. Order is preserved, because the
    5' and 3' partners define different transcripts. What this collapses is one
    directed pair reported repeatedly, once per alternative breakpoint or
    isoform. Each point here is a single caller, so the orientation it reports is
    its own claim; the figures that pool calls across callers merge the two
    orientations instead. Unparseable identifiers keep their raw string."""
    if not isinstance(fid, str) or not fid:
        return None
    parts = fid.replace('::', ':').split(':')
    return '::'.join(parts) if len(parts) >= 2 and all(parts) else fid


df['fusion_key'] = df['fusion.ens.gene.id'].map(_oriented)


def per_count(sub, name):
    """Distinct oriented fusions per (library, tool, replicate). Counting rows
    would confound how many fusions a library recovers with how verbosely each
    caller writes its output, which differs by up to 2.6-fold between libraries
    in this figure."""
    return (sub.dropna(subset=['fusion_key'])
              .groupby(KEYS)['fusion_key'].nunique().reset_index(name=name))


# One row per (library, tool, replicate) with total / known / novel counts.
merged = (per_count(df, 'total')
          .merge(per_count(df[df.Discovery.isin(KNOWN_TYPES)], 'known'), how='left')
          .merge(per_count(df[~df.Discovery.isin(KNOWN_TYPES)], 'novel'), how='left')
          .fillna({'known': 0.0, 'novel': 0.0}))

# (column, panel label, y-scale, ylim)
PANELS = [
    ('total', 'Total fusions',     'log',    (0.7, 2e6)),
    ('known', 'Known DB recovery', 'linear', (0,   40)),
    ('novel', 'Putative novel',    'log',    (0.7, 2e6)),
]


# ---------- Figure ----------
n_panels = len(PANELS)
n_libs = len(LIB_ORDER)

fig = plt.figure(figsize=(10.6, 8.6))
fig.patch.set_facecolor('white')

LEFT, RIGHT = 0.100, 0.975
TOP, BOT = 0.880, 0.085
GAP = 0.052
PANEL_H = (TOP - BOT - GAP * (n_panels - 1)) / n_panels


for pi, (col, label, scale, ylim) in enumerate(PANELS):
    y_pos = TOP - (pi + 1) * PANEL_H - pi * GAP
    ax = fig.add_axes([LEFT, y_pos, RIGHT - LEFT, PANEL_H])

    ax.set_yscale(scale)
    ax.set_ylim(*ylim)
    ax.set_xlim(-0.5, n_libs - 0.5)

    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(1.2)
        ax.spines[sp].set_color(NEUTRAL['ink'])
    ax.grid(True, axis='y', color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
    ax.tick_params(axis='y', labelsize=9, color=NEUTRAL['ink'])

    box_data = []
    for plat, lib, _l in LIB_ORDER:
        vals = merged[(merged.Platform == plat)
                      & (merged.library_type == lib)][col].values.astype(float)
        if scale == 'log':
            vals = np.where(vals == 0, 0.7, vals)
        box_data.append(vals)

    ax.boxplot(box_data, positions=range(n_libs), widths=0.55,
               whis=(0, 100), showfliers=False, patch_artist=True, zorder=3,
               medianprops=dict(color='black', linewidth=1.5),
               whiskerprops=dict(color=NEUTRAL['ink'], linewidth=0.9),
               capprops=dict(color=NEUTRAL['ink'], linewidth=0.9),
               boxprops=dict(facecolor='#F4F6F8', edgecolor=NEUTRAL['ink'],
                             linewidth=0.9))

    rng = np.random.default_rng(seed=pi * 17 + 3)
    for li, (plat, lib, _l) in enumerate(LIB_ORDER):
        sub = merged[(merged.Platform == plat) & (merged.library_type == lib)]
        for _, r in sub.iterrows():
            v = r[col]
            v = max(v, 0.7) if scale == 'log' else v
            ax.scatter(li + rng.uniform(-0.13, 0.13), v, s=46,
                       c=ALL_TOOL_COLORS.get(r['Algorithm'], NEUTRAL['muted']),
                       edgecolors='black', linewidths=0.4, alpha=0.9,
                       zorder=5, clip_on=False)

    ax.set_xticks(range(n_libs))
    if pi == n_panels - 1:
        ax.set_xticklabels([l for _, _, l in LIB_ORDER],
                           fontsize=9.2, color=NEUTRAL['body'])
    else:
        ax.set_xticklabels([])
    ax.tick_params(axis='x', length=0)

    unit = 'calls, log' if scale == 'log' else 'calls'
    ax.set_ylabel(f'{label}\n({unit})', fontsize=10.4, fontweight='bold',
                  color=NEUTRAL['ink'], labelpad=4)


# ---------- Title + tool legend ----------
fig.text(0.035, 0.955, 'Huh7 fusion-call distributions across libraries',
         ha='left', va='top', fontsize=12.8, fontweight='semibold',
         color=NEUTRAL['ink'])

TOOL_LEG = [t for t in ['CTAT-LR-Fusion', 'JAFFAL', 'Genion', 'FusionSeeker',
                        'GFSeeker', 'LongGF', 'JAFFA-direct', 'Arriba', 'STAR-Fusion']
            if t in set(df['Algorithm'])]
ax_leg = fig.add_axes([0, 0, 1, 1])
ax_leg.set_axis_off()
ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1)
leg_y = 0.912
leg_left, leg_right = 0.040, 0.965
spacing = (leg_right - leg_left) / len(TOOL_LEG)
for i, tool in enumerate(TOOL_LEG):
    lx = leg_left + i * spacing + 0.004
    ax_leg.scatter(lx, leg_y, s=80, c=ALL_TOOL_COLORS[tool],
                   edgecolors='white', linewidths=0.7,
                   transform=fig.transFigure, zorder=10, clip_on=False)
    ax_leg.text(lx + 0.011, leg_y, tool, ha='left', va='center',
                fontsize=8.4, color=NEUTRAL['body'], transform=fig.transFigure)


save_figure(fig, 'Fig3C', OUTDIR)
plt.close(fig)
