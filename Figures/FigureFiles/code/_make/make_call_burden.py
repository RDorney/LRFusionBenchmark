"""Generate call_burden.tsv: the number of DISTINCT fusions (gene pairs) detected
per Huh7 library for each caller, kept at spanning.reads >= 2 OR spanning.pairs >= 2
(the benchmark's read-support convention, matching the Huh7 figures).

A "fusion" is a distinct gene pair (fusion.ens.gene.id). Multiple call records for
the same pair differ mainly by breakpoint (breakpoint uncertainty / splicing) and
are collapsed to one; the representative record is the highest-read-support one.
n_fusions therefore counts distinct fusions, not call records.

Usage:
    python make_call_burden.py [output.tsv]     # default: call_burden.tsv beside this script
"""
import sys
from pathlib import Path
import pandas as pd
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, FUSIONS


DISCOVERY = FUSIONS / 'fusions_readsupport_Huh7_discovery.tsv.gz'
OUT = DATA / 'call_burden.tsv'

d = pd.read_csv(DISCOVERY, sep='\t',
               usecols=['Platform', 'library_type', 'RNA_sample', 'Algorithm',
                        'fusion.ens.gene.id', 'spanning.reads', 'spanning.pairs'],
               low_memory=False)
reads = pd.to_numeric(d['spanning.reads'], errors='coerce').fillna(0)
pairs = pd.to_numeric(d['spanning.pairs'], errors='coerce').fillna(0)
d = d[(reads >= 2) | (pairs >= 2)]

g = (d.groupby(['Algorithm', 'Platform', 'library_type', 'RNA_sample'])['fusion.ens.gene.id']
       .nunique().reset_index(name='n_fusions')
       .sort_values(['Algorithm', 'Platform', 'library_type', 'RNA_sample']))

g.to_csv(OUT, sep='\t', index=False)
print(f'wrote {OUT} ({len(g)} rows)')
