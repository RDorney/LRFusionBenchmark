"""Generate call_burden.tsv: the number of DISTINCT fusions (gene pairs) detected
per Huh7 library for each caller, kept at spanning.reads >= 2 OR spanning.pairs >= 2
(the benchmark's read-support convention, matching the Huh7 figures).

A fusion is a distinct ORIENTED gene pair. The 5' and 3' partners define different
transcripts, so A::B and B::A count separately; what is collapsed is one directed
pair reported repeatedly, once per alternative breakpoint or isoform. Degenerate
multi-segment calls such as B::A::B are excluded. Both rules come from the shared
loader so this table cannot drift from the figures that read it.

Usage:
    python make_call_burden.py [output.tsv]     # default: call_burden.tsv beside this script
"""
import sys
from pathlib import Path
import pandas as pd
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA
from _data import load, _fusion_key


OUT = DATA / 'call_burden.tsv'

# load() applies the read-support convention and drops degenerate calls.
d = load()
d = d.assign(fusion_key=d['fusion.ens.gene.id'].map(_fusion_key)).dropna(subset=['fusion_key'])

g = (d.groupby(['Algorithm', 'Platform', 'library_type', 'RNA_sample'])['fusion_key']
       .nunique().reset_index(name='n_fusions')
       .sort_values(['Algorithm', 'Platform', 'library_type', 'RNA_sample']))

g.to_csv(OUT, sep='\t', index=False)
print(f'wrote {OUT} ({len(g)} rows)')
