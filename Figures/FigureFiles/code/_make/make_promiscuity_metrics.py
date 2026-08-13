"""Build the per-library fusion-promiscuity metrics behind Figure 3F.

Each library (library_type x Platform x RNA_sample x Algorithm) is an undirected
gene network: nodes are genes, edges are distinct partnerships. A promiscuous gene
is a hub with many partners, which is a false-positive / chimeric-artefact signature.

Edge extraction, in full:
  * a record may name more than two genes (107 records name three or four), so it is
    expanded to ALL distinct gene pairs, C(n, 2);
  * self-fusions (A::A) are dropped, since a gene is not its own partner and a
    self-loop would otherwise be scored as two partners;
  * gene pairs are canonicalised (sorted) so A::B and B::A are one partnership, then
    deduplicated per library.

Read-support filter: kept fusions need `spanning.reads >= 2 OR spanning.pairs >= 2`,
which is the definition of a fusion call used throughout the benchmark. It is a no-op
for CTAT-LR-Fusion, FusionSeeker, GFSeeker, Genion and LongGF, trivial for JAFFAL, and
material only for JAFFA-direct.
"""
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, FUSIONS


DISC = FUSIONS / 'fusions_readsupport_Huh7_discovery.tsv.gz'
OUT = DATA
KEYS = ['library_type', 'Platform', 'RNA_sample', 'Algorithm']


def edges_of(x):
    """All distinct, non-self gene pairs named by one record."""
    if not isinstance(x, str):
        return []
    gs = [g for g in x.replace('::', ':').split(':') if g and g not in ('NA', 'nan')]  # drop unmapped (NA) genes
    return [tuple(sorted(p)) for p in combinations(gs, 2) if p[0] != p[1]]


def load(disc=DISC):
    d = pd.read_csv(disc, sep='\t',
                    usecols=KEYS + ['fusion.ens.gene.id', 'fusion.gene.id',
                                    'spanning.reads', 'spanning.pairs'],
                    low_memory=False)
    r = pd.to_numeric(d['spanning.reads'], errors='coerce').fillna(0)
    p = pd.to_numeric(d['spanning.pairs'], errors='coerce').fillna(0)
    d['keep'] = (r >= 2) | (p >= 2)
    # FusionSeeker's fusion.ens.gene.id is NA for ~all its rows: an upstream bug in
    # CheckKnownFusions.R sets it from the join-suffixed fusion.gene.id.x (which is NA)
    # rather than fusion.gene.id (which holds the real ENSG pair). Coalesce so FusionSeeker
    # is not silently dropped. No-op for every other tool, whose fusion.ens.gene.id is
    # populated.
    fid = d['fusion.ens.gene.id'].where(d['fusion.ens.gene.id'].notna(), d['fusion.gene.id'])
    d['edges'] = fid.map(edges_of)
    return d


def metrics(d, filtered):
    rows = []
    src = d[d.keep] if filtered else d
    for key, g in src.groupby(KEYS, sort=False):
        E = {e for lst in g['edges'] for e in lst}
        deg = Counter()
        for a, b in E:
            deg[a] += 1
            deg[b] += 1
        v = np.array(sorted(deg.values(), reverse=True)) if deg else np.array([0])
        total = v.sum() if v.sum() else 1
        lib, plat, rep, alg = key
        rows.append({
            'library_type': lib, 'Platform': plat, 'RNA_sample': rep, 'Algorithm': alg,
            'n_genes': float(len(deg)), 'n_pairs': float(len(E)),
            'mean_partners': round(float(v.mean()), 3),
            'median_partners': float(np.median(v)),
            'max_partners': float(v.max()),
            'prop_ge2': round(float((v >= 2).mean()), 4),
            'prop_ge3': round(float((v >= 3).mean()), 4),
            'prop_ge5': round(float((v >= 5).mean()), 4),
            'top1_share': round(float(v[:1].sum() / total), 4),
            'top5_share': round(float(v[:5].sum() / total), 4),
            'top10_share': round(float(v[:10].sum() / total), 4),
            'prep': f"{plat} {lib.replace('_', ' ')}",
        })
    return pd.DataFrame(rows).sort_values(KEYS).reset_index(drop=True)


def gene_degrees(d, filtered):
    rows = []
    src = d[d.keep] if filtered else d
    for key, g in src.groupby(KEYS, sort=False):
        E = {e for lst in g['edges'] for e in lst}
        deg = Counter()
        for a, b in E:
            deg[a] += 1
            deg[b] += 1
        lib, plat, rep, alg = key
        for gene, dg in deg.items():
            rows.append((lib, plat, rep, alg, gene, dg))
    return (pd.DataFrame(rows, columns=KEYS + ['Gene', 'degree'])
              .sort_values(KEYS + ['Gene'], kind='stable').reset_index(drop=True))


if __name__ == '__main__':
    disc = DISC
    if '--disc' in sys.argv:
        disc = Path(sys.argv[sys.argv.index('--disc') + 1])
        print(f'reading discovery table: {disc}')
    d = load(disc)

    if '--validate' in sys.argv:
        ref_path = OUT / sys.argv[sys.argv.index('--validate') + 1]
        ref = pd.read_csv(ref_path, sep='\t').sort_values(KEYS).reset_index(drop=True)
        got = metrics(d, filtered=False)
        cols = [c for c in ref.columns if c in got.columns]
        same = ref[cols].equals(got[cols])
        print(f'unfiltered rebuild vs {ref_path.name}: {len(got)} rows, identical = {same}')
        if not same:
            for c in cols:
                if not ref[c].equals(got[c]):
                    bad = (ref[c] != got[c])
                    print(f'  column {c}: {bad.sum()} rows differ')
                    print(pd.DataFrame({'ref': ref[c][bad], 'got': got[c][bad]}).head().to_string())
        sys.exit(0 if same else 1)

    m = metrics(d, filtered=True)
    m.to_csv(OUT / 'promiscuity_lib_metrics.tsv', sep='\t', index=False)
    print(f'wrote promiscuity_lib_metrics.tsv  ({len(m)} rows, reads>=2 | pairs>=2)')
    gd = gene_degrees(d, filtered=True)
    gd.to_csv(OUT / 'promiscuity_gene_degrees.tsv', sep='\t', index=False)
    print(f'wrote promiscuity_gene_degrees.tsv ({len(gd):,} rows)')
