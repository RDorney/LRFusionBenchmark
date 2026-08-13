"""Figure 2A. Overlap of recalled fusions across the six callers.

Each fusion recovered by at least one tool at 100Gb / 95% identity belongs to a
combination of tools. The dot matrix shows the membership of each combination,
the top bars show its size, and the left bars show each tool's total recall.

A fusion counts once, in the exact combination of tools that recovered it, and
the columns are ordered by bar height, ascending.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig2'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()


# ----- Intersections -----
ovl = pd.read_csv(DATA / 'fusions_in_n_algorithms.tsv', sep='\t')
sub = ovl[(ovl.depth == '100GB') & (ovl.Sequence_Identity == '95%') &
          (ovl.control == 'positive') & (ovl.recall_category == 'True_Recall')]

fusion_to_tools = sub.groupby('fusion.gene.id')['Algorithm'].apply(frozenset)

# One column per observed tool combination, ordered by size, ascending.
pattern_counts = fusion_to_tools.value_counts().sort_values(ascending=True)
patterns = list(pattern_counts.index)
heights = [int(pattern_counts[p]) for p in patterns]

set_sizes = {tool: int(sub[sub.Algorithm == tool]['fusion.gene.id'].nunique())
             for tool in TOOL_ORDER}


def render(y_label, title, basename):
    fig = plt.figure(figsize=(10.5, 5.6))
    fig.patch.set_facecolor('white')
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    VAL_RIGHT_x   = 0.035
    BAR_LEFT_x    = 0.050
    BAR_RIGHT_x   = 0.215
    TOOL_LBL_LEFT = 0.225
    MAIN_x0       = 0.345
    MAIN_x1       = 0.975
    TOP_y0        = 0.455
    TOP_y1        = 0.810
    DOT_y0        = 0.150
    DOT_y1        = 0.445

    n_cols = len(patterns)
    col_w = (MAIN_x1 - MAIN_x0) / n_cols
    bar_w_frac = 0.66

    TOP_MAX = int(np.ceil(max(heights) / 5) * 5)

    def y_top_frac(v):
        return TOP_y0 + v / TOP_MAX * (TOP_y1 - TOP_y0)

    # Top bars: count per combination
    INT_BAR_COLOR = '#3F6770'
    for j, n in enumerate(heights):
        x_mid = MAIN_x0 + (j + 0.5) * col_w
        bar_w = col_w * bar_w_frac
        x0 = x_mid - bar_w / 2
        y_top = y_top_frac(n)
        ax.add_patch(Rectangle((x0, TOP_y0), bar_w, y_top - TOP_y0,
                               transform=fig.transFigure,
                               facecolor=INT_BAR_COLOR, edgecolor='white',
                               linewidth=0.3, zorder=3))
        ax.text(x_mid, y_top + 0.006, str(n),
                transform=fig.transFigure,
                ha='center', va='bottom',
                fontsize=7.2, color=NEUTRAL['body'])

    top_ticks = list(range(0, TOP_MAX + 1, max(5, TOP_MAX // 5)))
    for tv in top_ticks:
        y = y_top_frac(tv)
        ax.plot([MAIN_x0, MAIN_x1], [y, y],
                transform=fig.transFigure,
                color=NEUTRAL['grid'], linewidth=0.4, zorder=1)
        ax.text(MAIN_x0 - 0.006, y, str(tv),
                transform=fig.transFigure,
                ha='right', va='center',
                fontsize=8.0, color=NEUTRAL['body'])
    ax.text(MAIN_x0 - 0.034, (TOP_y0 + TOP_y1) / 2, y_label,
            transform=fig.transFigure,
            rotation=90, ha='center', va='center',
            fontsize=9.0, color=NEUTRAL['body'])

    # Dot matrix
    n_tools = len(TOOL_ORDER)
    row_h = (DOT_y1 - DOT_y0) / n_tools
    DOT_OFF_COLOR = '#DCDCDC'
    LINE_COLOR    = '#4A4A4A'

    for r, tool in enumerate(TOOL_ORDER):
        y_bottom = DOT_y0 + (n_tools - 1 - r) * row_h
        if r % 2 == 0:
            ax.add_patch(Rectangle(
                (MAIN_x0, y_bottom), MAIN_x1 - MAIN_x0, row_h,
                transform=fig.transFigure,
                facecolor='#FAFAF6', edgecolor='none', zorder=0))

    for j, pat in enumerate(patterns):
        x_mid = MAIN_x0 + (j + 0.5) * col_w
        on_ys = []
        for r, tool in enumerate(TOOL_ORDER):
            y_mid = DOT_y0 + (n_tools - 1 - r) * row_h + row_h / 2
            is_on = tool in pat
            ax.scatter(
                x_mid, y_mid, s=64,
                c=TOOL_COLORS[tool] if is_on else DOT_OFF_COLOR,
                edgecolors='white', linewidths=0.8,
                transform=fig.transFigure, zorder=4, clip_on=False,
            )
            if is_on:
                on_ys.append(y_mid)
        if len(on_ys) >= 2:
            ax.plot([x_mid, x_mid], [min(on_ys), max(on_ys)],
                    transform=fig.transFigure,
                    color=LINE_COLOR, linewidth=2.0, zorder=3,
                    solid_capstyle='round')

    # Left set-size bars, growing left, with tool name between bar and matrix
    max_set = max(set_sizes.values())
    SET_MAX = int(np.ceil(max_set / 20) * 20)
    bar_range = BAR_RIGHT_x - BAR_LEFT_x

    def x_bar_tip(v):
        return BAR_RIGHT_x - v / SET_MAX * bar_range

    ax.text((BAR_LEFT_x + BAR_RIGHT_x) / 2, DOT_y1 + 0.012, 'Total recalled',
            transform=fig.transFigure,
            ha='center', va='bottom',
            fontsize=8.6, color=NEUTRAL['muted'])

    for r, tool in enumerate(TOOL_ORDER):
        y_bottom = DOT_y0 + (n_tools - 1 - r) * row_h
        y_mid = y_bottom + row_h / 2
        n = set_sizes[tool]
        bar_h = row_h * 0.55
        x_tip = x_bar_tip(n)
        ax.add_patch(Rectangle((x_tip, y_mid - bar_h / 2),
                               BAR_RIGHT_x - x_tip, bar_h,
                               transform=fig.transFigure,
                               facecolor=TOOL_COLORS[tool], edgecolor='white',
                               linewidth=0.3, zorder=3))
        ax.text(x_tip - 0.005, y_mid, str(n),
                transform=fig.transFigure,
                ha='right', va='center',
                fontsize=8.2, color=NEUTRAL['body'])
        ax.text(TOOL_LBL_LEFT, y_mid, tool,
                transform=fig.transFigure,
                ha='left', va='center',
                fontsize=8.8, color=NEUTRAL['body'])

    ax.text(0.030, 0.945, title,
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=12.5, fontweight='semibold', color=NEUTRAL['ink'])

    save_figure(fig, basename, OUTDIR)
    plt.close(fig)


render('Fusions detected by exactly this set',
       'Recalled-fusion overlap across callers',
       'Fig2A')
