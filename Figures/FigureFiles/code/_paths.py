"""Shared paths for the manuscript figure code.

Everything resolves relative to this file, so the package runs unchanged from any
location once the repository is cloned.

    from _paths import DATA, FUSIONS, OUT
    df = pd.read_csv(DATA / 'stat_summary.tsv', sep='\t')
"""

from pathlib import Path

# .../Figures/FigureFiles/code/_paths.py -> .../FigureFiles -> .../Figures -> repository root
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parent.parent

DATA = REPO / 'Figures' / 'Input_dataframes'   # every input table
FUSIONS = REPO / 'Fusions'                     # Huh7 discovery tables
OUT = ROOT / 'output'                          # rendered SVG and PNG

OUT.mkdir(parents=True, exist_ok=True)
