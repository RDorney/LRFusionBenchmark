"""Figure 1F. Breakpoint-position accuracy as a cumulative distribution.

Two panels (Breakpoint 1 and Breakpoint 2). Each curve is the ECDF of the
absolute distance between a caller's predicted breakpoint and the simulated
truth, on a log x-axis. Reference lines mark 100 bp, 1 kb and 1 Mb. At the
100 bp line each caller's percentage within tolerance is labelled, a star marks
callers that pass the TOST equivalence test, and a red triangle flags a worst
miss of at least 1 Mb. Conditions: 100Gb depth, 95% mean identity.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig1'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

set_style()

df = pd.read_csv(DATA / 'absolute_Breakpoint_Accuracy.tsv', sep='\t')
df = df[(df.seq_depth == '100GB') & (df.seq_id == '95%')
        & (df.Breakpoint_number.isin(['Breakpoint_1', 'Breakpoint_2']))].copy()
df['abs_d'] = df['Distance_from_Simulated_Breakpoint'].abs().astype(float)
df.loc[df.abs_d == 0, 'abs_d'] = 0.5

TOST_PASS = {('LongGF', 'Breakpoint_1'), ('Genion', 'Breakpoint_2')}
BPS    = ['Breakpoint_1', 'Breakpoint_2']
LABELS = {'Breakpoint_1': 'Breakpoint 1', 'Breakpoint_2': 'Breakpoint 2'}

REF_100 = 100
REF_1K  = 1_000
REF_1M  = 1_000_000


def fmt_bp_short(v):
    if v < 1000:   return f'{int(round(v))} bp'
    if v < 1e6:    return f'{v/1e3:.1f} kb'
    return f'{v/1e6:.1f} Mb'


# ---------- Figure ----------
fig = plt.figure(figsize=(11.2, 6.8))
fig.patch.set_facecolor('white')

L, R = 0.060, 0.985
B, T = 0.110, 0.745
GAP  = 0.075
W = (R - L - GAP) / 2

axes = []
for i, bp in enumerate(BPS):
    ax = fig.add_axes([L + i * (W + GAP), B, W, T - B])
    axes.append((ax, bp))


for ax, bp in axes:
    ax.set_xscale('log')
    ax.set_xlim(0.4, 5e8)
    ax.set_ylim(0, 1.05)

    # Reference verticals, with labels above the plot
    ax.axvline(REF_100, color='#1A1A1A', linewidth=1.8, zorder=2)
    ax.text(REF_100, 1.10, '±100 bp', ha='center', va='bottom',
            fontsize=10.6, fontweight='bold', color=NEUTRAL['ink'])

    ax.axvline(REF_1K, color='#9BB3C9', linewidth=0.8,
               linestyle=(0, (3, 3)), zorder=1)
    ax.text(REF_1K, 1.10, '1 kb', ha='center', va='bottom',
            fontsize=9.4, color='#5F7287')

    ax.axvline(REF_1M, color='#B2182B', linewidth=1.0,
               linestyle=(0, (3, 3)), zorder=1)
    ax.text(REF_1M, 1.10, '1 Mb', ha='center', va='bottom',
            fontsize=10.2, fontweight='bold', color='#B2182B')

    # Spines / grid
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(1.4)
        ax.spines[sp].set_color(NEUTRAL['ink'])
    ax.grid(True, color=NEUTRAL['grid'], linewidth=0.4, zorder=0)
    ax.tick_params(axis='both', labelsize=9.5, color=NEUTRAL['ink'])
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f'{int(round(v*100))}%'))

    # Per-tool ECDFs
    rows = []
    for tool in TOOL_ORDER:
        sub = df[(df.Algorithm == tool) & (df.Breakpoint_number == bp)]
        if not len(sub):
            continue
        x = np.sort(sub.abs_d.values)
        n = len(x)
        y = np.arange(1, n + 1) / n
        color = TOOL_COLORS[tool]

        ax.step(np.concatenate([[0.4], x]),
                np.concatenate([[0],   y]),
                where='post', color=color,
                linewidth=2.4, alpha=0.95, zorder=3)

        # % at 100bp
        pct = (sub.abs_d <= REF_100).sum() / n
        rows.append({'tool': tool, 'pct': pct, 'color': color,
                     'worst': x[-1], 'tost': (tool, bp) in TOST_PASS})

    # --- Dots / stars at 100 bp intersection (Q2 + Q4 combined) ---
    for r in rows:
        if r['tost']:
            ax.scatter([REF_100], [r['pct']], s=360, marker='*',
                       c=r['color'], edgecolors='black', linewidths=1.8,
                       zorder=8, clip_on=False)
        else:
            ax.scatter([REF_100], [r['pct']], s=85, c=r['color'],
                       edgecolors='black', linewidths=1.4,
                       zorder=7, clip_on=False)

    # --- Stack % labels with leader lines (avoid overlap) ---
    rows_sorted = sorted(rows, key=lambda r: -r['pct'])
    MIN_SEP = 0.07
    ys_disp = [r['pct'] for r in rows_sorted]
    for j in range(1, len(ys_disp)):
        if ys_disp[j - 1] - ys_disp[j] < MIN_SEP:
            ys_disp[j] = ys_disp[j - 1] - MIN_SEP
    # If bottommost was pushed below 0, shift the whole stack up uniformly
    if ys_disp[-1] < 0:
        bump = -ys_disp[-1]
        ys_disp = [y + bump for y in ys_disp]

    label_x = REF_100 * 2.6
    for r, y_disp in zip(rows_sorted, ys_disp):
        ax.annotate(f'{int(round(r["pct"]*100))}%',
                    xy=(REF_100, r['pct']),
                    xytext=(label_x, y_disp),
                    textcoords='data',
                    ha='left', va='center',
                    fontsize=10.0, fontweight='bold',
                    color=r['color'], zorder=9,
                    arrowprops=dict(arrowstyle='-', color=r['color'],
                                    alpha=0.55, lw=0.8,
                                    shrinkA=4, shrinkB=4,
                                    connectionstyle='arc3,rad=0'),
                    bbox=dict(boxstyle='round,pad=0.20',
                              facecolor='white', edgecolor='none',
                              alpha=0.92))

    # --- Catastrophe markers (Q3) ---
    catastrophes = [r for r in rows if r['worst'] >= REF_1M]
    catastrophes.sort(key=lambda r: r['worst'])
    for k, r in enumerate(catastrophes):
        worst = r['worst']
        # Triangle at curve endpoint (Y=1.0)
        ax.scatter([worst], [1.0], s=170, marker='v',
                   c='#B2182B', edgecolors='black', linewidths=1.0,
                   zorder=8, clip_on=False)
        # Stagger labels when several coincide. Text in the tool colour, red
        # border for the catastrophe alert.
        y_off = -10 - k * 14
        ax.annotate(f'worst miss: {fmt_bp_short(worst)}',
                    xy=(worst, 1.0),
                    xytext=(-8, y_off), textcoords='offset points',
                    ha='right', va='top',
                    fontsize=9.0, fontweight='bold',
                    color=r['color'], zorder=9,
                    bbox=dict(boxstyle='round,pad=0.20',
                              facecolor='white',
                              edgecolor='#B2182B', linewidth=1.1,
                              alpha=0.95))

    ax.set_title(LABELS[bp], fontsize=12.4, fontweight='semibold',
                 color=NEUTRAL['ink'], loc='left', pad=46)
    ax.set_xlabel('|Distance from simulated breakpoint|  (bp, log)',
                  fontsize=10.8, fontweight='bold',
                  color=NEUTRAL['ink'], labelpad=4)

axes[0][0].set_ylabel('Cumulative fraction of fusion predictions',
                      fontsize=11, fontweight='bold', color=NEUTRAL['ink'])
axes[1][0].tick_params(axis='y', labelleft=False)


# ---------- Figure title + key ----------
ax_lab = fig.add_axes([0, 0, 1, 1])
ax_lab.set_axis_off()
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)

ax_lab.text(0.035, 0.962,
            'Breakpoint-position accuracy (ECDF)',
            transform=fig.transFigure,
            ha='left', va='top',
            fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'])

# Tool colour legend (scatter dots, evenly spaced)
leg_y = 0.890
leg_left, leg_right = 0.040, 0.730
spacing = (leg_right - leg_left) / len(TOOL_ORDER)
for i, tool in enumerate(TOOL_ORDER):
    lx = leg_left + i * spacing + 0.005
    ax_lab.scatter(lx, leg_y, s=90, c=TOOL_COLORS[tool],
                   edgecolors='white', linewidths=0.7,
                   transform=fig.transFigure, zorder=4, clip_on=False)
    ax_lab.text(lx + 0.012, leg_y, tool,
                transform=fig.transFigure,
                ha='left', va='center',
                fontsize=8.7, color=NEUTRAL['body'])

# Marker key, right of tool legend
key_y = leg_y
ax_lab.scatter(0.760, key_y, s=240, marker='*', c='#7B7B7B',
               edgecolors='black', linewidths=1.4,
               transform=fig.transFigure, zorder=4, clip_on=False)
ax_lab.text(0.772, key_y, 'TOST pass',
            transform=fig.transFigure,
            ha='left', va='center', fontsize=8.6, color=NEUTRAL['body'])
ax_lab.scatter(0.870, key_y, s=120, marker='v', c='#B2182B',
               edgecolors='black', linewidths=0.9,
               transform=fig.transFigure, zorder=4, clip_on=False)
ax_lab.text(0.882, key_y, 'worst ≥ 1 Mb',
            transform=fig.transFigure,
            ha='left', va='center', fontsize=8.6, color=NEUTRAL['body'])

save_figure(fig, 'Fig1F', OUTDIR)
plt.close(fig)
