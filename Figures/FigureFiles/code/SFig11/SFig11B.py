"""Supplementary Figure 11B. Sense-antisense fusion calls in Huh7.

One panel per library prep (PCR cDNA, direct cDNA, direct RNA); within each, a
boxplot per platform of the number of sense-antisense calls across callers and
replicates, with the individual caller values overlaid as coloured points on a
log count axis. Only CTAT-LR-Fusion, FusionSeeker, GFSeeker, JAFFA-direct and
JAFFAL produce these calls; Genion and LongGF make none. Calls are kept at a
minimum of two spanning reads or read pairs.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import FUSIONS, OUT
OUTDIR = OUT / 'SFig11'
from palette import TOOL_COLORS, SHORT_READ_COLORS, NEUTRAL, set_style, save_figure

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

sa = df[df.fusionType == 'Sense-Antisense']
agg = (sa.groupby(['library_type', 'Platform', 'RNA_sample', 'Algorithm'])['fusion_key']
         .nunique().reset_index(name='count'))

LIBRARIES   = ['PCR_cDNA', 'direct_cDNA', 'direct_RNA']
LIB_TITLES  = {'PCR_cDNA': 'PCR cDNA',
               'direct_cDNA': 'Direct cDNA',
               'direct_RNA':  'Direct RNA'}
PLATFORMS   = ['Illumina', 'PacBio', 'ONT']
ALL_TOOL_COLORS = {**TOOL_COLORS, **SHORT_READ_COLORS}
TOOL_ORDER_LOCAL = [t for t in ['CTAT-LR-Fusion', 'FusionSeeker', 'GFSeeker', 'JAFFAL',
                                'JAFFA-direct', 'Arriba', 'STAR-Fusion']
                    if t in set(agg['Algorithm'])]


# ---------- Figure ----------
fig = plt.figure(figsize=(11.6, 5.2))
fig.patch.set_facecolor('white')

L, R = 0.060, 0.985
B, T = 0.150, 0.745
GAP  = 0.038
W = (R - L - GAP * (len(LIBRARIES) - 1)) / len(LIBRARIES)

axes = []
for i, lib in enumerate(LIBRARIES):
    ax = fig.add_axes([L + i * (W + GAP), B, W, T - B])
    axes.append((ax, lib))


for ax, lib in axes:
    ax.set_yscale('log')
    ax.set_ylim(0.7, 1.5e4)
    ax.set_xlim(-0.5, len(PLATFORMS) - 0.5)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(1.2)
        ax.spines[sp].set_color(NEUTRAL['ink'])
    ax.grid(True, axis='y', color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
    ax.tick_params(axis='y', labelsize=9.5, color=NEUTRAL['ink'])

    sub = agg[agg.library_type == lib]
    # Collect distributions per Platform for boxplot
    box_data = []
    for plat in PLATFORMS:
        vals = sub[sub.Platform == plat]['count'].values
        box_data.append(vals if len(vals) else np.array([np.nan]))

    # Box: whiskers extend to data min/max (matches paper convention)
    bp = ax.boxplot(box_data, positions=range(len(PLATFORMS)),
                    widths=0.55, whis=(0, 100), showfliers=False,
                    patch_artist=True, zorder=3,
                    medianprops=dict(color='black', linewidth=1.6),
                    whiskerprops=dict(color=NEUTRAL['ink'], linewidth=0.9),
                    capprops=dict(color=NEUTRAL['ink'], linewidth=0.9),
                    boxprops=dict(facecolor='#F4F6F8',
                                  edgecolor=NEUTRAL['ink'],
                                  linewidth=0.9))

    # Overlay points coloured by tool
    rng = np.random.default_rng(seed=42)
    for plat_i, plat in enumerate(PLATFORMS):
        plat_sub = sub[sub.Platform == plat]
        for _, r in plat_sub.iterrows():
            jit = rng.uniform(-0.10, 0.10)
            ax.scatter(plat_i + jit, r['count'], s=58,
                       c=ALL_TOOL_COLORS[r['Algorithm']],
                       edgecolors='black', linewidths=0.5,
                       alpha=0.92, zorder=5, clip_on=False)

    # Set x-tick labels AFTER boxplot (boxplot would overwrite them)
    ax.set_xticks(range(len(PLATFORMS)))
    ax.set_xticklabels(PLATFORMS, fontsize=9.5, color=NEUTRAL['body'])

    ax.set_title(LIB_TITLES[lib], fontsize=11.4, fontweight='semibold',
                 color=NEUTRAL['ink'], loc='left', pad=8)

axes[0][0].set_ylabel('Sense-antisense calls  (count, log)',
                      fontsize=10.6, fontweight='bold',
                      color=NEUTRAL['ink'], labelpad=4)
for ax, _ in axes[1:]:
    ax.tick_params(axis='y', labelleft=False)


# ---------- Title block + legend ----------
ax_lab = fig.add_axes([0, 0, 1, 1])
ax_lab.set_axis_off()
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)

ax_lab.text(0.035, 0.952,
            'Sense-antisense fusion calls in Huh7',
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=12.6, fontweight='semibold', color=NEUTRAL['ink'])

leg_y = 0.892
leg_left, leg_right = 0.040, 0.940
spacing = (leg_right - leg_left) / len(TOOL_ORDER_LOCAL)
for i, tool in enumerate(TOOL_ORDER_LOCAL):
    lx = leg_left + i * spacing + 0.005
    ax_lab.scatter(lx, leg_y, s=90, c=ALL_TOOL_COLORS[tool],
                   edgecolors='white', linewidths=0.7,
                   transform=fig.transFigure, zorder=4, clip_on=False)
    ax_lab.text(lx + 0.012, leg_y, tool,
                transform=fig.transFigure,
                ha='left', va='center',
                fontsize=8.7, color=NEUTRAL['body'])


save_figure(fig, 'SFig11B', OUTDIR)
plt.close(fig)
