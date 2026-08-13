"""Figure 4D. Huh7 known fusion-call similarity across library replicates."""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _similarity_helper import render_similarity
from _paths import OUT

render_similarity(fusion_filter='known', basename='Fig4D',
                  title='Huh7 known fusion-call similarity across libraries',
                  outdir=OUT / 'Fig4')
