"""Supplementary Figure 13. Per-read alignment identity by Huh7 library.

One box per replicate, grouped into four library preps and sharing a single
identity axis so the platform hierarchy is directly comparable. Identity is the
per-read match fraction (aligned_length - NM) / aligned_length from the genome
BAM, with introns (N) and clips excluded. Boxes are coloured by platform and
whiskers span the 1st-99th percentile. Values are a 150k-read sample per library.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'SFig13'
from palette import NEUTRAL, set_style, save_figure

set_style()

df = pd.read_csv(DATA / 'huh7_read_identity_sample.tsv.gz', sep='\t')

# Library order and titles match the read-length figure so the two panels stack
# consistently.
SAMPLES = [
    ('ONT',    'direct_RNA'),
    ('ONT',    'direct_cDNA'),
    ('ONT',    'PCR_cDNA'),
    ('PacBio', 'PCR_cDNA'),
]
SAMPLE_TITLES = {
    ('ONT',    'direct_RNA'):  'ONT\nDirect RNA',
    ('ONT',    'direct_cDNA'): 'ONT\nDirect cDNA',
    ('ONT',    'PCR_cDNA'):    'ONT\nPCR cDNA',
    ('PacBio', 'PCR_cDNA'):    'PacBio\nPCR cDNA',
}
PLATFORM_COLORS = {'ONT': '#1F6BA0', 'PacBio': '#C73E5C'}
REPS = ['B1', 'B2']

# Assemble the eight boxes and their x positions: two replicates per cluster,
# with a wider gap between clusters.
PAIR_GAP    = 0.90
CLUSTER_GAP = 1.55
BOX_W       = 0.66

data, positions, colors, rep_labels = [], [], [], []
cluster_centres, cluster_titles = [], []
x = 0.0
for platform, lib in SAMPLES:
    xs = []
    for rep in REPS:
        vals = df[(df.Platform == platform)
                  & (df.library_type == lib)
                  & (df.replicate == rep)]['identity'].values
        data.append(vals)
        positions.append(x)
        colors.append(PLATFORM_COLORS[platform])
        rep_labels.append(rep)
        xs.append(x)
        x += PAIR_GAP
    cluster_centres.append(np.mean(xs))
    cluster_titles.append(SAMPLE_TITLES[(platform, lib)])
    x += CLUSTER_GAP - PAIR_GAP

medians = [float(np.median(v)) for v in data]
q95s    = [float(np.percentile(v, 95)) for v in data]

# Whiskers span the 5th-95th percentile (central 90% of reads). A few preps have
# a long low-identity tail below this that would otherwise flatten the axis.
WHIS = (5, 95)
p5_min = min(float(np.percentile(v, 5)) for v in data)
y_lo = np.floor((p5_min - 1.5) / 2) * 2
y_hi = 100.9


# ---------- Figure ----------
fig = plt.figure(figsize=(13.4, 5.6))
fig.patch.set_facecolor('white')

L, R = 0.070, 0.985
B, T = 0.150, 0.755
ax = fig.add_axes([L, B, R - L, T - B])

bp = ax.boxplot(
    data, positions=positions, widths=BOX_W,
    patch_artist=True, whis=WHIS, showfliers=False,
    medianprops=dict(color='white', linewidth=1.8),
    whiskerprops=dict(linewidth=1.2),
    capprops=dict(linewidth=1.2),
    zorder=3,
)
for patch, box, whisk, cap, color in zip(
        bp['boxes'], bp['boxes'],
        [bp['whiskers'][i:i + 2] for i in range(0, len(bp['whiskers']), 2)],
        [bp['caps'][i:i + 2] for i in range(0, len(bp['caps']), 2)],
        colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.78)
    patch.set_edgecolor(color)
    patch.set_linewidth(1.2)
    for w in whisk:
        w.set_color(color)
    for c in cap:
        c.set_color(color)

# Median value just above each box's upper whisker
for pos, med, q95, color in zip(positions, medians, q95s, colors):
    ax.text(pos, q95 + 0.10, f'{med:.1f}',
            ha='center', va='bottom', fontsize=7.6, fontweight='bold',
            color=color, alpha=0.95)

# B1 / B2 tick labels under each box
ax.set_xticks(positions)
ax.set_xticklabels(rep_labels, fontsize=8.6, color=NEUTRAL['body'])
ax.set_xlim(positions[0] - BOX_W, positions[-1] + BOX_W)

# Cluster (library prep) labels below the replicate labels
for centre, title in zip(cluster_centres, cluster_titles):
    ax.text(centre, -0.085, title, transform=ax.get_xaxis_transform(),
            ha='center', va='top', fontsize=9.6, fontweight='semibold',
            color=NEUTRAL['ink'], linespacing=1.1)

ax.set_ylim(y_lo, y_hi)
ax.set_ylabel('Read alignment identity  (%)',
              fontsize=10.8, fontweight='bold', color=NEUTRAL['ink'], labelpad=4)

for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
for sp in ('left', 'bottom'):
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color(NEUTRAL['ink'])
ax.grid(True, axis='y', color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
ax.tick_params(axis='y', labelsize=9, color=NEUTRAL['ink'])
ax.tick_params(axis='x', length=0)


# ---------- Title block + legend ----------
ax_lab = fig.add_axes([0, 0, 1, 1])
ax_lab.set_axis_off()
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)

ax_lab.text(0.035, 0.952,
            'Huh7 per-read alignment identity',
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])

leg_y = 0.862
ax_lab.scatter(0.040, leg_y, s=90, c=PLATFORM_COLORS['ONT'],
               edgecolors='white', linewidths=0.7,
               transform=fig.transFigure, zorder=4, clip_on=False)
ax_lab.text(0.052, leg_y, 'ONT',
            transform=fig.transFigure,
            ha='left', va='center', fontsize=9.0, color=NEUTRAL['body'])
ax_lab.scatter(0.120, leg_y, s=90, c=PLATFORM_COLORS['PacBio'],
               edgecolors='white', linewidths=0.7,
               transform=fig.transFigure, zorder=4, clip_on=False)
ax_lab.text(0.132, leg_y, 'PacBio',
            transform=fig.transFigure,
            ha='left', va='center', fontsize=9.0, color=NEUTRAL['body'])


save_figure(fig, 'SFig13', OUTDIR)
plt.close(fig)
