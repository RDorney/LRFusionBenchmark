"""Figure 3E and Supplementary Figure 12. Fusion-partner circos plots on the hg38 cytoband ideogram, one panel per caller and library.

Generates one Circos per long-read caller x library preparation (Huh7, replicate
B1 only), 6 callers x 4 preps = 24 figures. The four GFSeeker panels are the
main-text series; the full set is the supplementary material.

Every panel is rendered identically so they are directly comparable: chromosomes
are the real hg38 Giemsa-banded ideogram, each arc is a distinct fusion gene pair
(spanning reads >= 2 or spanning pairs >= 2), and each panel colours its own top 5
promiscuous genes by partner count (ties broken alphabetically). A gene needs at
least HUB_MIN partners to qualify, matching the >= 2 definition of promiscuity used
throughout the analysis: a gene with a single partner is one fusion, not a hub. The
legend reports each gene's partner count, so a clean caller is self-evident: its top
gene has a handful of partners rather than the hundreds seen in promiscuous panels,
and a panel with no legend has no promiscuous gene at all.
"""
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import networkx as nx
from pycirclize import Circos
from pycirclize.utils import load_eukaryote_example_dataset

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from palette import NEUTRAL, save_figure
from _paths import DATA, FUSIONS, OUT
OUTDIR = OUT / 'Fig3'

DISC = FUSIONS / 'fusions_readsupport_Huh7_discovery.tsv.gz'
SPANS = DATA / 'ref_annot.gtf.gene_spans'
GRID_BASE = OUT

REP = 'B1'
N_HUBS = 5                      # top 5 among the qualifying genes
HUB_MIN = 2                     # >= 2 partners to count as promiscuous (analysis-wide definition)
HUB_PALETTE = ['#E67C6E', '#F0A94E', '#6FA3CE', '#8FC77F', '#B48FC7']
WIN = 1_500_000                 # +/- window (bp) giving each link ribbon width


def link_style(n):
    """Background and hub link styling, scaled to how many links a panel draws.
    The direct preps draw only a handful of fusions, which disappear at the faint
    hairline weight the dense PCR-cDNA panels need, so sparse panels use darker,
    thicker, more opaque links. Returns (bg_colour, bg_alpha, bg_lw, hub_alpha,
    hub_lw); the densest bucket keeps the original faint styling unchanged."""
    if n <= 60:
        return '#565E68', 0.70, 0.55, 0.85, 0.65
    if n <= 250:
        return '#6E767E', 0.42, 0.38, 0.80, 0.45
    if n <= 800:
        return '#868D95', 0.20, 0.20, 0.78, 0.28
    return '#9DA4AB', 0.11, 0.10, 0.75, 0.18

# (platform, library_type, short name for the file, pretty prep label)
# Two modes: the default long-read grid (6 callers x 4 ONT/PacBio preps) and the
# Illumina short-read set (3 callers x Illumina PCR-cDNA), written to Grid/Illumina.
# Usage: python N10_circos_grid.py [longread|illumina] [caller]
MODE = sys.argv[1] if len(sys.argv) > 1 else 'longread'
if MODE == 'illumina':
    CALLERS = ['STAR-Fusion', 'Arriba', 'JAFFA-direct']   # quickest first; JAFFA-direct is the 507k-arc monster
    LIBS = [('Illumina', 'PCR_cDNA', 'Illumina_PCRcDNA', 'Illumina PCR-cDNA')]
    OUTDIR = GRID_BASE / 'Illumina'
else:
    CALLERS = ['CTAT-LR-Fusion', 'JAFFAL', 'Genion', 'FusionSeeker', 'GFSeeker', 'LongGF']
    LIBS = [
        ('ONT', 'direct_RNA', 'ONT_directRNA', 'ONT direct RNA'),
        ('ONT', 'direct_cDNA', 'ONT_directcDNA', 'ONT direct cDNA'),
        ('ONT', 'PCR_cDNA', 'ONT_PCRcDNA', 'ONT PCR-cDNA'),
        ('PacBio', 'PCR_cDNA', 'PacBio_PCRcDNA', 'PacBio PCR-cDNA'),
    ]
    OUTDIR = GRID_BASE

ONLY = sys.argv[2] if len(sys.argv) > 2 else None   # optional: render one caller only
if ONLY:
    CALLERS = [c for c in CALLERS if c == ONLY]

OUTDIR.mkdir(parents=True, exist_ok=True)

# ----- gene loci + symbols -----
gs = pd.read_csv(SPANS, sep='\t', header=None, usecols=[0, 1, 2, 3, 5],
                 names=['ensg_v', 'chrom', 'start', 'end', 'symbol'])
gs['ensg'] = gs['ensg_v'].str.split('.').str[0]
gs['mid'] = ((gs['start'] + gs['end']) / 2).astype(int)
POS = {r.ensg: (r.chrom, r.mid) for r in gs.itertuples()}
SYM = dict(zip(gs['ensg'], gs['symbol']))

# ----- discovery table, read once -----
d = pd.read_csv(DISC, sep='\t',
                usecols=['Platform', 'library_type', 'RNA_sample', 'Algorithm',
                         'fusion.ens.gene.id', 'spanning.reads', 'spanning.pairs'],
                low_memory=False)
d['r'] = pd.to_numeric(d['spanning.reads'], errors='coerce').fillna(0)
d['p'] = pd.to_numeric(d['spanning.pairs'], errors='coerce').fillna(0)
d = d[(d.r >= 2) | (d.p >= 2)]


def edges_of(x):
    """All distinct gene pairs named by one record. Most records name two genes, but
    107 in this table name three or four, so keeping only the first two would silently
    drop partnerships. Self-fusions (A::A) are not partnerships and are dropped. This
    matches how promiscuity_lib_metrics.tsv was built."""
    if not isinstance(x, str):
        return []
    gs = [g for g in x.replace('::', ':').split(':') if g and g not in ('NA', 'nan')]  # drop unmapped (NA) genes
    return [tuple(sorted(p)) for p in combinations(gs, 2) if p[0] != p[1]]


d['edges'] = d['fusion.ens.gene.id'].map(edges_of)
d = d[d.RNA_sample == REP]

chr_bed_file = DATA / 'hg38_chr.bed'
cytoband_file = DATA / 'hg38_cytoband.tsv'
if not (chr_bed_file.exists() and cytoband_file.exists()):
    chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")


def hub_labels(top):
    """Legend labels 'SYMBOL (partners)', with the locus added where a symbol is
    ambiguous. Pseudoautosomal genes (CD99, ASMTL, PLCXD1 ...) carry separate chrX
    and chrY gene IDs in GENCODE v44, and 330 symbols in the annotation map to more
    than one gene ID, so the symbol alone can name two distinct nodes."""
    syms = [SYM.get(g, g) for g, _ in top]
    out = []
    for (g, dg), s in zip(top, syms):
        if syms.count(s) > 1 and g in POS:
            s = f'{s} {POS[g][0]}'
        out.append(f'{s} ({dg})')
    return out


def make_figure(caller, plat, lib, short, prep_label):
    circos = Circos.initialize_from_bed(chr_bed_file, space=3)
    circos.add_cytoband_tracks((95, 100), cytoband_file)
    CHR_SIZE = {s.name: s.size for s in circos.sectors}
    for sector in circos.sectors:
        sector.text(sector.name.replace('chr', ''), r=110, size=9, color=NEUTRAL['body'])

    def region(ensg):
        if ensg not in POS:
            return None
        chrom, mid = POS[ensg]
        if chrom not in CHR_SIZE:
            return None
        return (chrom, max(0, mid - WIN), min(int(CHR_SIZE[chrom]), mid + WIN))

    sub = d[(d.Algorithm == caller) & (d.Platform == plat) & (d.library_type == lib)]
    edges = sorted({e for lst in sub['edges'] for e in lst})
    edges = [e for e in edges if region(e[0]) is not None and region(e[1]) is not None]
    G = nx.Graph()
    G.add_edges_from(edges)
    deg = dict(G.degree())
    ranked = sorted(deg.items(), key=lambda kv: (-kv[1], SYM.get(kv[0], kv[0])))
    top = [(g, dg) for g, dg in ranked if dg >= HUB_MIN][:N_HUBS]
    HUB_BY_ENSG = {g: HUB_PALETTE[i] for i, (g, _) in enumerate(top)}
    max_deg = max(deg.values()) if deg else 0

    bg_col, bg_alpha, bg_lw, hub_alpha, hub_lw = link_style(len(edges))
    hub_links = []
    for a, b in edges:
        ra, rb = region(a), region(b)
        if ra is None or rb is None:
            continue
        col = HUB_BY_ENSG.get(a) or HUB_BY_ENSG.get(b)
        if col:
            hub_links.append((ra, rb, col))
        else:
            circos.link(ra, rb, color=bg_col, alpha=bg_alpha, lw=bg_lw)
    for ra, rb, col in hub_links:
        circos.link(ra, rb, color=col, alpha=hub_alpha, lw=hub_lw)

    fig = circos.plotfig(dpi=100)
    fig.suptitle(f'{caller}  ·  Huh7 {prep_label} ({REP})',
                 x=0.5, y=1.03, fontsize=12, fontweight='semibold', color=NEUTRAL['ink'])
    fig.text(0.5, 0.995,
             f'{G.number_of_nodes():,} genes  ·  {G.number_of_edges():,} fusions  ·  '
             f'max {max_deg} partners',
             ha='center', va='top', fontsize=9.5, color=NEUTRAL['muted'])
    if top:
        h = [Line2D([0], [0], color=HUB_PALETTE[i], lw=3, label=lab)
             for i, lab in enumerate(hub_labels(top))]
        fig.legend(handles=h, loc='lower center', ncol=len(top), frameon=False,
                   fontsize=9, bbox_to_anchor=(0.5, -0.02))
    panel = {'JAFFAL': 'Fig3E', 'CTAT-LR-Fusion': 'SFig12A', 'FusionSeeker': 'SFig12B',
             'Genion': 'SFig12C', 'GFSeeker': 'SFig12D', 'LongGF': 'SFig12E'}[caller]
    name = f'{panel}_{short}'
    save_figure(fig, name, GRID_BASE / ('Fig3' if panel == 'Fig3E' else 'SFig12'))
    plt.close(fig)
    print(f'{name}: {G.number_of_nodes()} genes, {G.number_of_edges()} fusions, '
          f'max {max_deg}, hubs {[(SYM.get(g, g), dg) for g, dg in top]}', flush=True)


n = 0
for caller in CALLERS:
    for plat, lib, short, label in LIBS:
        make_figure(caller, plat, lib, short, label)
        n += 1
        print(f'--- {n}/{len(CALLERS) * len(LIBS)} done ---', flush=True)
print('done')
