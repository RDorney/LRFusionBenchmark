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
KNOWN_TYPES = ['Known', 'Reverse Known', 'contains Known 1.2', 'contains Known 2.3']


def load(min_read_support=2):
    """Calls kept at spanning.reads >= n (long read) or spanning.pairs >= n (Illumina)."""
    df = pd.read_csv(DATA, sep='\t', low_memory=False)
    reads = pd.to_numeric(df['spanning.reads'], errors='coerce').fillna(0)
    pairs = pd.to_numeric(df['spanning.pairs'], errors='coerce').fillna(0)
    df['spanning.reads'] = pd.to_numeric(df['spanning.reads'], errors='coerce')
    df = df[(reads >= min_read_support) | (pairs >= min_read_support)]
    return df


def _canon(fid):
    """Canonical gene-pair id: sort A:B so A::B and B::A collapse to one."""
    if not isinstance(fid, str):
        return None
    parts = fid.replace('::', ':').split(':')
    if len(parts) < 2:
        return None
    return ':'.join(sorted(parts))


def _canon_set(series):
    """Apply canonical normalisation and return the unique set."""
    s = set()
    for v in series.dropna():
        c = _canon(v)
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
