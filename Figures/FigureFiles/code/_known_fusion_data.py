"""
Data layer for the Huh7 known-fusion recovery figure (Figure 4C).

Fusions are matched to caller output by Ensembl gene-ID pair rather than by
symbol string. Gene IDs are stable across symbol renames and independent of gene
order, so a fusion is credited to a library whichever symbol version and
orientation the caller reports. Column labels use the modern (v110) symbols.

Reversed-orientation duplicates (A::B and B::A share one gene-ID key) are
collapsed to a single fusion, and by default only fusions detected in at least
one library are returned, so the figure shows fusions that genuinely appear in
the data rather than the full curated list.

Source list: KnownHuh7Fusion_ah.csv, curated and Ensembl v110 annotated.
"""

import numpy as np
import pandas as pd

from _paths import REPO, FUSIONS

AH_CSV    = REPO / 'Known_Fusions' / 'KnownHuh7Fusion_ah.csv'
DISCOVERY = FUSIONS / 'fusions_readsupport_Huh7_discovery.tsv.gz'

# Pipeline flags marking a call as a known fusion.
KNOWN_TYPES = ['Known', 'Reverse Known', 'contains Known 1.2', 'contains Known 2.3']

# Upset rows (library types), top to bottom.
LIBRARY_ROWS = [
    ('ONT',      'direct_RNA',  'ONT direct RNA'),
    ('ONT',      'direct_cDNA', 'ONT direct cDNA'),
    ('ONT',      'PCR_cDNA',    'ONT PCR cDNA'),
    ('PacBio',   'PCR_cDNA',    'PacBio PCR cDNA'),
    ('Illumina', 'PCR_cDNA',    'Illumina PCR cDNA'),
]

# Literature sample origin mapped to display category.
_CATEGORY = {'HEPG2': 'HCC', 'Fibrolamellar-HCC': 'HCC', 'HCC': 'HCC', 'HUH7': 'HUH7'}


def _clean(v):
    return v.strip() if isinstance(v, str) else v


def _modern(vname, orig):
    """Modern v110 symbol, with fallback to the original symbol when the v110
    lookup is empty or NA. NaN is truthy, so the test is explicit."""
    v = _clean(vname)
    if isinstance(v, str) and v:
        return v
    return _clean(orig)


def _parse_ens_pair(s):
    """Parse fusion.ens.gene.id into the frozenset of its two Ensembl gene IDs.

    Order and version agnostic: ENSG_A::ENSG_B equals ENSG_B::ENSG_A, and a
    trailing version suffix is stripped. Returns None when fewer than two
    distinct ENSG tokens are present.
    """
    if not isinstance(s, str):
        return None
    ids = set()
    for tok in s.replace('::', ':').split(':'):
        tok = tok.strip().split('.')[0]
        if tok.startswith('ENSG'):
            ids.add(tok)
    return frozenset(ids) if len(ids) >= 2 else None


def load_known_fusions():
    """Curated literature fusions with a modern label and an Ensembl-ID key.

    Columns: label (modern A::B), orig_name (as written in the literature CSV),
    category (HCC or HUH7), key (frozenset of the two Ensembl gene IDs, or None
    when unavailable).
    """
    lit = pd.read_csv(AH_CSV)
    recs = []
    for _, r in lit.iterrows():
        gx, gy = _clean(r['GENEID.x']), _clean(r['GENEID.y'])
        mx = _modern(r['v110_name.x'], r['Gene1'])
        my = _modern(r['v110_name.y'], r['Gene2'])
        valid = (isinstance(gx, str) and gx.startswith('ENSG')
                 and isinstance(gy, str) and gy.startswith('ENSG') and gx != gy)
        recs.append({
            'label':     f'{mx}::{my}',
            'orig_name': _clean(r['Known_Gene_Fusion']),
            'category':  _CATEGORY.get(_clean(r['Sample']), _clean(r['Sample'])),
            'key':       frozenset({gx, gy}) if valid else None,
        })
    return pd.DataFrame(recs)


def _dedup_by_key(kf):
    """Collapse rows that share a gene-ID key (reversed orientation or synonym),
    keeping the first occurrence. Rows without a key are left untouched."""
    seen, rows = set(), []
    for _, r in kf.iterrows():
        k = r['key']
        if k is not None:
            if k in seen:
                continue
            seen.add(k)
        rows.append(r)
    return pd.DataFrame(rows).reset_index(drop=True)


def _det_map(frame):
    """Map (Platform, library_type) to the set of detected gene-ID frozensets."""
    return {k: set(g['_key'])
            for k, g in frame.groupby(['Platform', 'library_type'], sort=False)}


def _presence(kf, detmap):
    M = np.zeros((len(LIBRARY_ROWS), len(kf)), dtype=bool)
    keys = kf['key'].tolist()
    for j, key in enumerate(keys):
        if key is None:
            continue
        for i, (plat, lib, _lbl) in enumerate(LIBRARY_ROWS):
            M[i, j] = key in detmap.get((plat, lib), set())
    return M


def prepare(read_support=None, restrict_known=False, detected_only=True,
            verbose=False):
    """Build the data for the known-fusion renderer.

    Parameters
    ----------
    read_support : int or None
        Keep only calls with spanning.reads >= n or spanning.pairs >= n.
        None applies no support filter.
    restrict_known : bool
        If True, count only calls the pipeline flagged as known (Discovery in
        KNOWN_TYPES). If False, match on gene ID against the curated list.
    detected_only : bool
        If True (default), drop fusions detected by no library so only fusions
        present in the data are shown.

    Returns
    -------
    all_fusions : list of str
        Modern labels, HCC block then HUH7 block, each ordered by detection
        count then name.
    category_idx : list of str
        'HCC' or 'HUH7' per column, aligned with all_fusions.
    M : numpy.ndarray of bool
        Libraries (rows) by fusions (columns) presence matrix.
    """
    kf = _dedup_by_key(load_known_fusions())

    usecols = ['Platform', 'library_type', 'fusion.ens.gene.id',
               'Discovery', 'spanning.reads', 'spanning.pairs']
    disc = pd.read_csv(DISCOVERY, sep='\t', usecols=usecols, low_memory=False)

    if read_support:
        reads = pd.to_numeric(disc['spanning.reads'], errors='coerce').fillna(0)
        pairs = pd.to_numeric(disc['spanning.pairs'], errors='coerce').fillna(0)
        disc = disc[(reads >= read_support) | (pairs >= read_support)]

    disc = disc.assign(_key=disc['fusion.ens.gene.id'].map(_parse_ens_pair))
    disc = disc[disc['_key'].notna()]

    detect_all   = _det_map(disc)
    detect_known = _det_map(disc[disc['Discovery'].isin(KNOWN_TYPES)])
    detect = detect_known if restrict_known else detect_all

    M_full = _presence(kf, detect)
    n_det = M_full.sum(axis=0)

    labels = kf['label'].tolist()
    cats = kf['category'].tolist()
    keep = [j for j in range(len(kf)) if (n_det[j] >= 1 or not detected_only)]

    order = []
    for cat in ['HCC', 'HUH7']:
        idx = [j for j in keep if cats[j] == cat]
        idx.sort(key=lambda j: (-int(n_det[j]), str(labels[j])))
        order.extend(idx)

    all_fusions  = [labels[j] for j in order]
    category_idx = [cats[j] for j in order]
    M = M_full[:, order]

    if verbose:
        _report(kf, n_det, detected_only)

    return all_fusions, category_idx, M


def _report(kf, n_det, detected_only):
    """Print an audit of the label mapping and detection counts."""
    print(f'[known-fusion data]  source = {AH_CSV.name}   '
          f'({len(kf)} fusions after dedup)')
    print(f'  categories: {kf["category"].value_counts().to_dict()}')
    n_shown = int((n_det >= 1).sum())
    print(f'  detected by >=1 library: {n_shown}   '
          f'(detected_only={detected_only})')
    nd = [kf.iloc[j]['label'] for j in range(len(kf)) if n_det[j] == 0]
    print(f'  {len(nd)} fusion(s) detected by no library'
          f'{" (dropped)" if detected_only else ""}')


if __name__ == '__main__':
    fusions, cats, M = prepare(verbose=True)
    print(f'  -> {len(fusions)} columns drawn, matrix {M.shape}')
