"""Supplementary Figure 10. Huh7 fusion calls by fusion type across libraries.

Six stacked panels, one per fusion type. In each panel the x-axis is the five
Huh7 libraries and the box shows the spread across tools and replicates of the
number of DISTINCT fusions of that type, with individual tool values overlaid as
dots. Panels use independent y-scales because per-type counts span several
orders of magnitude. Calls are filtered to a minimum read support of two.

A fusion is an ORIENTED Ensembl gene-ID pair. The 5' and 3' partners define
different transcripts with different reading frames, so A::B and B::A are
counted separately; what is collapsed is the repeated reporting of one directed
pair, which a caller emits once per alternative breakpoint or isoform. Counting
rows instead would confound how many fusions a library recovers with how
verbosely each caller writes its output, and that differs by up to fivefold
within a single panel (CTAT-LR-Fusion reports 3,682 sense-antisense rows for 575
distinct pairs in ONT direct cDNA).

The figures that pool calls across callers (Figure 4C to 4E, Supplementary
Figures 9 and 11B) merge the two orientations instead, because callers disagree
on which partner they place first for the same underlying event. That does not
apply here: each point is a single caller, so the orientation it reports is its
own claim and is preserved. Within one caller, library and replicate the two
choices are almost identical anyway, differing by 0 to 58 pairs for every
long-read caller. The exception is JAFFA-direct, which reports both orientations
for 230,664 gene pairs in the Illumina libraries; those reciprocal calls are
shown rather than merged away.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import OUT
OUTDIR = OUT / 'SFig10'
from palette import TOOL_COLORS, SHORT_READ_COLORS, NEUTRAL, set_style, save_figure
from _data import load, LIB_ORDER

set_style()

df = load()  # minimum read support of two


def _oriented(fid):
    """Oriented gene-pair key, 5' partner first. Order is deliberately preserved.
    75 Arriba rows carry a half-empty id such as 'ENSG0000...::' because one
    partner did not map to an Ensembl gene, and no fallback column is populated
    for them. Those keep the raw string as their key rather than being dropped,
    so the call still counts once and identical strings still merge."""
    if not isinstance(fid, str) or not fid:
        return None
    parts = fid.replace('::', ':').split(':')
    if len(parts) == 2 and all(parts):
        return '::'.join(parts)
    return fid


df['fusion_key'] = df['fusion.ens.gene.id'].map(_oriented)

KEYS = ['Platform', 'library_type', 'Algorithm', 'RNA_sample']
# Every tool x library x replicate that produced at least one call. Per-type
# counts are zero-filled over this set so a tool that calls none of a type shows
# a genuine zero rather than dropping out.
base = df[KEYS].drop_duplicates()

# The Illumina short-read callers (JAFFA-direct, Arriba, STAR-Fusion) sit outside the
# long-read palette and carry their own tints. Without them they fell back to grey,
# which reads as LongGF (also grey) in the Illumina column.
ALL_TOOL_COLORS = {**TOOL_COLORS, **SHORT_READ_COLORS}

# (fusionType value, panel label, y-scale, ylim), ordered by call volume.
# SAGe is included: cis-splicing between adjacent genes is a recognised
# transcript class rather than an obvious artefact, it is discussed in the
# accompanying text, and it appears in no other figure. The obvious-artefact
# classes (Mitochondrial:Genomic, Mitochondrial:Mitochondrial,
# Self-Misalignment) remain in Supplementary Figure 11, which is dedicated to
# them, and tetra-fusion is too sparse to plot.
PANELS = [
    ('inter-chromosomal', 'Inter-chromosomal', 'log',    (0.7, 2e6)),
    ('intra-chromosomal', 'Intra-chromosomal', 'log',    (0.7, 1e5)),
    ('Sense-Antisense',   'Sense-antisense',   'log',    (0.7, 2e3)),
    ('read-through',      'Read-through',      'log',    (0.7, 1e3)),
    ('SAGe',              'SAGe',              'log',    (0.7, 2e2)),
    ('tri-fusion',        'Tri-fusion',        'linear', (0,   28)),
]


def type_counts(ftype):
    """Per (library, tool, replicate) count of distinct oriented fusions of one
    type, zero-filled over every tool x library x replicate that produced any
    call."""
    cnt = (df[df.fusionType == ftype]
             .dropna(subset=['fusion_key'])
             .groupby(KEYS)['fusion_key'].nunique()
             .reset_index(name='count'))
    return base.merge(cnt, how='left').fillna({'count': 0.0})


# ---------- Figure ----------
n_panels = len(PANELS)
n_libs = len(LIB_ORDER)

# The layout is written in inches and converted to figure fractions, so adding
# a panel grows the canvas instead of squeezing the existing panels or
# stretching the title block. At five panels these values reproduce the previous
# 10.6 x 12.6 in layout exactly.
PANEL_IN, GAP_IN, HEAD_IN, FOOT_IN = 1.75, 0.48, 1.11, 0.82
FIG_W = 10.6
FIG_H = HEAD_IN + n_panels * PANEL_IN + (n_panels - 1) * GAP_IN + FOOT_IN

fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor('white')

LEFT, RIGHT = 0.100, 0.975
TOP, BOT = 1 - HEAD_IN / FIG_H, FOOT_IN / FIG_H
GAP = GAP_IN / FIG_H
PANEL_H = PANEL_IN / FIG_H


for pi, (ftype, label, scale, ylim) in enumerate(PANELS):
    y_pos = TOP - (pi + 1) * PANEL_H - pi * GAP
    ax = fig.add_axes([LEFT, y_pos, RIGHT - LEFT, PANEL_H])
    m = type_counts(ftype)

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

    # Per-library counts across (tool, replicate)
    box_data = []
    for plat, lib, _l in LIB_ORDER:
        vals = m[(m.Platform == plat) & (m.library_type == lib)]['count'].values.astype(float)
        if scale == 'log':
            vals = np.where(vals == 0, 0.7, vals)   # tiny floor so zeros plot
        box_data.append(vals)

    ax.boxplot(box_data, positions=range(n_libs), widths=0.55,
               whis=(0, 100), showfliers=False, patch_artist=True, zorder=3,
               medianprops=dict(color='black', linewidth=1.5),
               whiskerprops=dict(color=NEUTRAL['ink'], linewidth=0.9),
               capprops=dict(color=NEUTRAL['ink'], linewidth=0.9),
               boxprops=dict(facecolor='#F4F6F8', edgecolor=NEUTRAL['ink'],
                             linewidth=0.9))

    # Overlay per-tool dots
    rng = np.random.default_rng(seed=pi * 17 + 3)
    for li, (plat, lib, _l) in enumerate(LIB_ORDER):
        sub = m[(m.Platform == plat) & (m.library_type == lib)]
        for _, r in sub.iterrows():
            v = r['count']
            v = max(v, 0.7) if scale == 'log' else v
            ax.scatter(li + rng.uniform(-0.13, 0.13), v, s=42,
                       c=ALL_TOOL_COLORS.get(r['Algorithm'], NEUTRAL['muted']),
                       edgecolors='black', linewidths=0.4, alpha=0.9,
                       zorder=5, clip_on=False)

    # Library labels on the bottom panel only
    ax.set_xticks(range(n_libs))
    if pi == n_panels - 1:
        ax.set_xticklabels([l for _, _, l in LIB_ORDER],
                           fontsize=9.2, color=NEUTRAL['body'])
    else:
        ax.set_xticklabels([])
    ax.tick_params(axis='x', length=0)

    unit = 'fusions, log' if scale == 'log' else 'fusions'
    ax.set_ylabel(f'{label}\n({unit})', fontsize=10.2, fontweight='bold',
                  color=NEUTRAL['ink'], labelpad=4)


# ---------- Title + tool legend ----------
fig.text(0.035, 1 - 0.315 / FIG_H,
         'Distinct Huh7 fusions by fusion type across libraries',
         ha='left', va='top', fontsize=12.8, fontweight='semibold',
         color=NEUTRAL['ink'])

TOOL_LEG = [t for t in ['CTAT-LR-Fusion', 'JAFFAL', 'Genion', 'FusionSeeker',
                        'GFSeeker', 'LongGF', 'JAFFA-direct', 'Arriba', 'STAR-Fusion']
            if t in set(df['Algorithm'])]
ax_leg = fig.add_axes([0, 0, 1, 1])
ax_leg.set_axis_off()
ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1)
leg_y = 1 - 0.693 / FIG_H
leg_left, leg_right = 0.040, 0.965
spacing = (leg_right - leg_left) / len(TOOL_LEG)
for i, tool in enumerate(TOOL_LEG):
    lx = leg_left + i * spacing + 0.004
    ax_leg.scatter(lx, leg_y, s=80, c=ALL_TOOL_COLORS[tool],
                   edgecolors='white', linewidths=0.7,
                   transform=fig.transFigure, zorder=10, clip_on=False)
    ax_leg.text(lx + 0.011, leg_y, tool, ha='left', va='center',
                fontsize=8.4, color=NEUTRAL['body'], transform=fig.transFigure)


save_figure(fig, 'SFig10', OUTDIR)
plt.close(fig)
