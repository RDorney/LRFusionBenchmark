"""Shared data loader for the Huh7 fusion-call figures."""

import pandas as pd

from _paths import FUSIONS

DATA = FUSIONS / 'fusions_readsupport_Huh7_discovery.tsv.gz'

LIB_ORDER = [
    ('ONT',      'direct_RNA',  'ONT\ndRNA'),
    ('ONT',      'direct_cDNA', 'ONT\ndcDNA'),
    ('ONT',      'PCR_cDNA',    'ONT\nPCR cDNA'),
    ('PacBio',   'PCR_cDNA',    'PacBio\nPCR cDNA'),
    ('Illumina', 'PCR_cDNA',    'Illumina\nPCR cDNA'),
]
# A call counts as a known-fusion recovery only when the gene order matches the
# literature. 'Reverse Known' has the same pair with the 5' and 3' partners
# swapped, which is a different transcript with a different reading frame, so it
# is deliberately absent here and falls on the novel side of the split. The
# 'contains Known' labels only ever applied to degenerate B::A::B calls, which
# are now excluded upstream, so they are gone too.
KNOWN_TYPES = ['Known']


def load(min_read_support=2):
    """Calls kept at spanning.reads >= n (long read) or spanning.pairs >= n (Illumina)."""
    df = pd.read_csv(DATA, sep='\t', low_memory=False)
    reads = pd.to_numeric(df['spanning.reads'], errors='coerce').fillna(0)
    pairs = pd.to_numeric(df['spanning.pairs'], errors='coerce').fillna(0)
    df['spanning.reads'] = pd.to_numeric(df['spanning.reads'], errors='coerce')
    df = df[(reads >= min_read_support) | (pairs >= min_read_support)]
    # Degenerate multi-segment calls such as B::A::B are excluded outright.
    df = df[~df['fusion.ens.gene.id'].map(is_degenerate)]
    return df


def _fusion_key(fid):
    """Identity of a fusion: its Ensembl gene IDs in the order reported.

    Orientation is PRESERVED. The 5' and 3' partners define different
    transcripts with different reading frames, so A::B and B::A are two
    different fusions and are never merged. A key does collapse one directed
    pair reported repeatedly, once per alternative breakpoint or isoform.
    Version suffixes are stripped. An identifier that does not yield two
    Ensembl IDs, such as the half-empty '::ENSG00000084674' that Arriba and
    JAFFA-direct emit when one partner fails to map, keeps its raw string so
    the call is keyed rather than dropped.
    """
    if not isinstance(fid, str) or not fid:
        return None
    ids = [t.split('.')[0] for t in fid.replace('::', ':').split(':')
           if t.startswith('ENSG')]
    return '::'.join(ids) if len(ids) >= 2 else fid


def is_degenerate(fid):
    """True for a multi-segment call that repeats a gene, such as B::A::B.

    That is not a genuine multi-gene fusion, it is one gene pair written three
    times over, so it is excluded from every figure rather than folded into
    the two-gene pair. A real three-gene fusion A::B::C keeps three distinct
    IDs and is retained.
    """
    if not isinstance(fid, str):
        return False
    ids = [t.split('.')[0] for t in fid.replace('::', ':').split(':')
           if t.startswith('ENSG')]
    return len(ids) >= 3 and len(set(ids)) < len(ids)


def _canon_set(series):
    """Unique fusion keys in a slice, orientation preserved."""
    s = set()
    for v in series.dropna():
        c = _fusion_key(v)
        if c is not None:
            s.add(c)
    return s


def per_library_counts(df):
    """Return one row per library × stage with unique-fusion counts.
    Uses canonical gene-pair IDs (A::B and B::A merged) to avoid
    double-counting fusions that different tools report in different orders.
    """
    rows = []
    for plat, lib, label in LIB_ORDER:
        sub = df[(df.Platform == plat) & (df.library_type == lib)]
        sub_k   = sub[sub.Discovery.isin(KNOWN_TYPES)]
        sub_rt  = sub_k[sub_k.Validation.fillna('').str.contains('RT-PCR')]
        sub_h   = sub_k[sub_k.Sample.fillna('') == 'HUH7']
        sub_hcc = sub_k[sub_k.Sample.fillna('').isin(
            ['HCC', 'HEPG2', 'Fibrolamellar-HCC'])]
        total_ids = _canon_set(sub['fusion.ens.gene.id'])
        known_ids = _canon_set(sub_k['fusion.ens.gene.id'])
        rt_ids    = _canon_set(sub_rt['fusion.ens.gene.id'])
        huh_ids   = _canon_set(sub_h['fusion.ens.gene.id'])
        hcc_ids   = _canon_set(sub_hcc['fusion.ens.gene.id'])
        hcc_only  = len(hcc_ids - huh_ids)
        rows.append({
            'platform': plat, 'library_type': lib, 'label': label,
            'total':   len(total_ids),
            'known':   len(known_ids),
            'rtpcr':   len(rt_ids),
            'huh7_lit': len(huh_ids),
            'hcc_only': hcc_only,
        })
    return pd.DataFrame(rows)
