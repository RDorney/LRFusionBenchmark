"""
Shared rendering for the fusion-call similarity heatmaps.

Each axis is one library replicate (platform, prep and replicate). Its fusion set
is the union across all callers, deduplicated to canonical gene-pair IDs so that
A::B and B::A collapse to one junction. Cells are the Jaccard similarity between
two replicates' sets, and rows and columns are reordered by average-linkage
clustering on 1 - Jaccard. Calls are kept at spanning.reads >= 2 or
spanning.pairs >= 2.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform

from palette import F1_CMAP, NEUTRAL, set_style, save_figure
from _paths import FUSIONS, OUT

DATA = FUSIONS / 'fusions_readsupport_Huh7_discovery.tsv.gz'
OUTDIR = OUT
KNOWN_TYPES = ['Known', 'Reverse Known', 'contains Known 1.2', 'contains Known 2.3']

# Short, readable prep names for the axis labels.
_LIB_SHORT = {'direct_RNA': 'dRNA', 'direct_cDNA': 'dcDNA', 'PCR_cDNA': 'PCR-cDNA'}


def _library_label(row):
    """Platform, short prep and replicate, e.g. 'ONT dRNA B1'. The discovery
    file has no combined library label, so it is built here."""
    lib = _LIB_SHORT.get(row['library_type'], row['library_type'])
    return f"{row['Platform']} {lib} {row['RNA_sample']}"


def _canon(fid):
    """Canonical gene-pair id: sort A:B so A::B and B::A collapse to one."""
    if not isinstance(fid, str):
        return None
    parts = fid.replace('::', ':').split(':')
    if len(parts) < 2:
        return None
    return ':'.join(sorted(parts))


def render_similarity(*, fusion_filter, basename, title, outdir=OUTDIR):
    """Render one Jaccard-similarity heatmap.

    fusion_filter in {'all', 'known', 'novel'} selects which discovery rows are
    used before building the per-replicate fusion sets.
    """
    set_style()
    df = pd.read_csv(DATA, sep='\t', low_memory=False)
    df['spanning.reads'] = pd.to_numeric(df['spanning.reads'], errors='coerce').fillna(0)
    df['spanning.pairs'] = pd.to_numeric(df['spanning.pairs'], errors='coerce').fillna(0)
    df = df[(df['spanning.reads'] >= 2) | (df['spanning.pairs'] >= 2)]

    if fusion_filter == 'known':
        df = df[df.Discovery.isin(KNOWN_TYPES)]
    elif fusion_filter == 'novel':
        df = df[df.Discovery == 'Putative Novel']
    elif fusion_filter != 'all':
        raise ValueError(f"unknown fusion_filter {fusion_filter!r}")

    df = df.assign(full_label=df.apply(_library_label, axis=1))

    sets = {}
    for lab, sub in df.groupby('full_label'):
        ids = set()
        for raw in sub['fusion.ens.gene.id'].dropna():
            c = _canon(raw)
            if c is not None:
                ids.add(c)
        sets[lab] = ids
    labels = sorted(sets.keys())
    n = len(labels)
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            a, b = sets[labels[i]], sets[labels[j]]
            u = len(a | b)
            J[i, j] = len(a & b) / u if u else 0

    # Hierarchical clustering on 1 - J
    D = 1 - J
    np.fill_diagonal(D, 0)
    Z = linkage(squareform(D, checks=False), method='average')
    order = leaves_list(Z)
    labels_ord = [labels[i] for i in order]
    J_ord = J[np.ix_(order, order)]

    # Heatmap kept square in display inches.
    FIG_W, FIG_H = 9.6, 7.2
    HM_INCHES = 4.6
    HM_W = HM_INCHES / FIG_W
    HM_H = HM_INCHES / FIG_H

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor('white')

    L = 0.300
    HM_B = 0.150
    DD_H = 0.095
    CBAR_GAP = 0.020
    CBAR_W = 0.018

    ax_top = fig.add_axes([L, HM_B + HM_H + 0.008, HM_W, DD_H])
    ax_hm = fig.add_axes([L, HM_B, HM_W, HM_H])
    ax_cb = fig.add_axes([L + HM_W + CBAR_GAP, HM_B, CBAR_W, HM_H])

    dendrogram(Z, ax=ax_top, orientation='top',
               color_threshold=0, above_threshold_color=NEUTRAL['rule'],
               no_labels=True)
    ax_top.set_xlim(0, n * 10)
    for sp in ax_top.spines.values():
        sp.set_visible(False)
    ax_top.tick_params(left=False, right=False, top=False, bottom=False,
                       labelleft=False, labelright=False,
                       labeltop=False, labelbottom=False)
    ax_top.set_xticks([]); ax_top.set_yticks([])

    im = ax_hm.imshow(J_ord, cmap=F1_CMAP, vmin=0, vmax=1,
                      aspect='equal', interpolation='nearest', origin='upper')
    ax_hm.set_xticks(range(n))
    ax_hm.set_yticks(range(n))
    ax_hm.set_xticklabels(labels_ord, rotation=38, ha='right',
                          fontsize=9.0, color=NEUTRAL['body'])
    ax_hm.set_yticklabels(labels_ord, fontsize=9.0, color=NEUTRAL['body'])

    for i in range(n):
        for j in range(n):
            v = J_ord[i, j]
            color = 'white' if v > 0.55 else NEUTRAL['ink']
            ax_hm.text(j, i, f'{v:.2f}', ha='center', va='center',
                       fontsize=8.4, color=color)

    ax_hm.set_xticks(np.arange(n + 1) - 0.5, minor=True)
    ax_hm.set_yticks(np.arange(n + 1) - 0.5, minor=True)
    ax_hm.grid(which='minor', color='white', linewidth=1.0)
    ax_hm.tick_params(which='minor', length=0)
    ax_hm.tick_params(which='major', length=0)
    for sp in ax_hm.spines.values():
        sp.set_visible(False)

    cb = fig.colorbar(im, cax=ax_cb, orientation='vertical')
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8.6, color=NEUTRAL['ink'], length=3)
    cb.set_label('Jaccard similarity', fontsize=9.6, fontweight='bold',
                 color=NEUTRAL['ink'], labelpad=8, rotation=270, va='bottom')

    fig.text(0.035, 0.955, title, ha='left', va='top',
             fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])

    save_figure(fig, basename, outdir)
    plt.close(fig)
