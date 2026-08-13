"""
Shared style and palette for the long-read fusion caller benchmark figures.

Design principles
-----------------
- Colour-blind safe and print safe.
- Sequential colour maps avoid rainbow; the categorical tool palette balances
  warm and cool hues.
- A single typography stack across all figures.
- Vector SVG and 600 dpi PNG export.

Usage
-----
    from palette import TOOL_ORDER, TOOL_COLORS, F1_CMAP, set_style
    set_style()
"""

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# ----- Tool ordering and categorical palette -----
# Order follows the manuscript ranking by F1 at 100GB / 95% identity:
# CTAT-LR-Fusion > JAFFAL > Genion > FusionSeeker > GFSeeker > LongGF.
TOOL_ORDER = [
    'CTAT-LR-Fusion',
    'JAFFAL',
    'Genion',
    'FusionSeeker',
    'GFSeeker',
    'LongGF',
]

# Categorical hues. Warm red for the top performer, cool and neutral for the
# rest. Chosen to stay distinct under deuteranopia and protanopia.
TOOL_COLORS = {
    'CTAT-LR-Fusion': '#D6604D',
    'JAFFAL':         '#4393C3',
    'Genion':         '#5AAE61',
    'FusionSeeker':   '#E08B2B',
    'GFSeeker':       '#9970AB',
    'LongGF':         '#7F8C8D',
}

# The Illumina short-read comparator (JAFFA-direct) sits outside the six
# long-read tools and uses this accent where it appears alongside them.
JAFFA_DIRECT_COLOR = '#7DBFD2'

# All three Illumina short-read comparators, coloured as a family distinct from the
# six long-read tools (JAFFA-direct light blue, Arriba magenta, STAR-Fusion teal).
SHORT_READ_COLORS = {
    'JAFFA-direct': '#7DBFD2',
    'Arriba':       '#CC6699',
    'STAR-Fusion':  '#3F9E9E',
}


# ----- Library-prep ordering + platform colours -----
# PREP_ORDER is the canonical top-to-bottom / left-to-right ordering of the
# real-data library preparations. Colour on the real-data figures is carried by
# PLATFORM, matching the read-length and per-read identity figures; this is the single
# platform palette for the paper.
PREP_ORDER = [
    'ONT direct RNA',
    'ONT direct cDNA',
    'ONT PCR-cDNA',
    'PacBio PCR-cDNA',
    'Illumina PCR-cDNA',
]
PLATFORM_COLORS = {
    'ONT':    '#1F6BA0',
    'PacBio': '#C73E5C',
}


def prep_key(platform, library_type):
    """Normalise a (platform, library_type) pair to a canonical PREP_ORDER key.

    Tolerates the spelling variants used across the input tables
    ('direct_RNA', 'direct RNA', 'PCR_cDNA', 'PCR-cDNA', ...).
    """
    plat = str(platform).strip()
    lib = str(library_type).replace('_', ' ').replace('-', ' ').strip().lower()
    if 'direct' in lib and 'rna' in lib:
        chem = 'direct RNA'
    elif 'direct' in lib and 'cdna' in lib:
        chem = 'direct cDNA'
    elif 'pcr' in lib:
        chem = 'PCR-cDNA'
    else:
        chem = str(library_type)
    return f'{plat} {chem}'


# ----- Sequential colour map for F1 / Precision / Recall heatmaps -----
# Cream to teal to deep navy, with monotonic luminance.
_F1_STOPS = [
    (0.00, '#F6EFE2'),
    (0.18, '#D8E0CE'),
    (0.36, '#9CC0B1'),
    (0.55, '#4F947F'),
    (0.78, '#2A6973'),
    (1.00, '#173E5C'),
]
F1_CMAP = LinearSegmentedColormap.from_list(
    'crest_warm', [c for _, c in _F1_STOPS], N=256
)

# Warmer ramp for absolute counts such as false positives.
_COUNT_STOPS = [
    (0.00, '#FBF3EC'),
    (0.30, '#F2C396'),
    (0.60, '#E07A47'),
    (1.00, '#7A2D1A'),
]
COUNT_CMAP = LinearSegmentedColormap.from_list(
    'ember', [c for _, c in _COUNT_STOPS], N=256
)

# Per-metric colour maps for the Precision / Recall / F1 stack. Each metric uses
# a different hue over a shared cream baseline. Blue, amber and teal stay
# distinct under deuteranopia and protanopia.
PRECISION_CMAP = LinearSegmentedColormap.from_list(
    'precision_blue',
    ['#F6EFE2', '#CBD9E5', '#7AA3C2', '#3973A0', '#0F3158'],
    N=256,
)
RECALL_CMAP = LinearSegmentedColormap.from_list(
    'recall_amber',
    ['#F6EFE2', '#EBD3B0', '#D8A267', '#A56826', '#5A3007'],
    N=256,
)
# F1 reuses F1_CMAP defined above.


# ----- Typography and figure defaults -----
FONT_STACK = ['Liberation Sans', 'Arial', 'DejaVu Sans']  # Arial-metric fallback

NEUTRAL = {
    'ink':      '#1A1A1A',   # title and strong text
    'body':     '#2E2E2E',   # body text
    'muted':    '#6B6B6B',   # captions and subtitles
    'rule':     '#C9C9C9',   # spines and fine rules
    'grid':     '#E8E8E8',   # gridlines
    'na_fill':  '#F0F0F0',   # NA cell background
    'na_text':  '#A5A5A5',   # NA dash
}


def set_style():
    """Apply the shared rcParams. Call once at the top of each figure script."""
    mpl.rcParams.update({
        # Fonts
        'font.family': 'sans-serif',
        'font.sans-serif': FONT_STACK,
        'font.size': 9.5,
        'axes.titlesize': 11.5,
        'axes.titleweight': 'semibold',
        'axes.labelsize': 9.5,
        'xtick.labelsize': 9.0,
        'ytick.labelsize': 9.0,
        'legend.fontsize': 8.5,
        'figure.titlesize': 12,

        # Colour of text and axes
        'text.color':       NEUTRAL['body'],
        'axes.labelcolor':  NEUTRAL['body'],
        'axes.edgecolor':   NEUTRAL['rule'],
        'xtick.color':      NEUTRAL['body'],
        'ytick.color':      NEUTRAL['body'],

        # Keep only the outer chrome
        'axes.spines.top':    False,
        'axes.spines.right':  False,
        'axes.linewidth':     0.6,
        'xtick.major.width':  0.6,
        'ytick.major.width':  0.6,
        'xtick.major.size':   3,
        'ytick.major.size':   3,

        # White background for print
        'figure.facecolor': 'white',
        'axes.facecolor':   'white',
        'savefig.facecolor': 'white',
        'savefig.edgecolor': 'none',
        'savefig.bbox':     'tight',
        'savefig.pad_inches': 0.08,

        # Vector and raster export defaults
        'pdf.fonttype':  42,       # TrueType, so text stays editable in a vector editor
        'ps.fonttype':   42,
        'svg.fonttype':  'none',   # keep text as <text> elements in the SVG
        'savefig.dpi':   600,
        'figure.dpi':    150,      # display dpi while iterating

        # Regular (non-italic) digits in math text
        'mathtext.default': 'regular',
    })


def auto_text_color(value, vmin, vmax, threshold=0.55):
    """Pick black or white text on a coloured cell from a luminance proxy."""
    if value is None:
        return NEUTRAL['na_text']
    t = (value - vmin) / max(vmax - vmin, 1e-9)
    return 'white' if t >= threshold else NEUTRAL['ink']


def save_figure(fig, basename, outdir):
    """Save an editable vector SVG and a 600 dpi PNG preview."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f'{basename}.svg')
    fig.savefig(outdir / f'{basename}.png', dpi=600)
    print(f'  saved {outdir / basename}.svg')
    print(f'  saved {outdir / basename}.png (600 dpi)')
