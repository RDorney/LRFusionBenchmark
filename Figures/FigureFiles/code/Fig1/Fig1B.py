"""Figure 1B. Computational performance across speed, memory and accuracy.

Speed, Memory and Accuracy, each with the six callers placed by RANK (best at
top). A tool-coloured line links every caller across the three axes, so the
crossing pattern shows the argument directly: a different tool leads each axis,
and no caller stays near the top of all three.

  Speed    rank by ascending user CPU minutes at 10 Gb (LongGF fastest).
  Memory   rank by ascending peak RAM at 10 Gb (Genion lightest).
  Accuracy rank by descending F1 at 10 Gb (CTAT-LR-Fusion most accurate).
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import DATA, OUT
OUTDIR = OUT / 'Fig1'
from palette import TOOL_ORDER, TOOL_COLORS, NEUTRAL, set_style, save_figure

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

set_style()

PERF = DATA / 'perf_value.tsv'


def fmt_time(m):
    return f"{m:.0f} min" if m < 90 else f"{m/60:.1f} h"


# ----- Load and rank -----
df = pd.read_csv(PERF, sep='\t').set_index('tool')
N = len(df)

# Best -> worst. Lower time / lower RAM / higher F1 is better.
speed_rank = df['user_min_10gb'].rank(method='min').astype(int)      # 1 = fastest
mem_rank   = df['ram_gb_10gb'].rank(method='min').astype(int)        # 1 = lightest
acc_rank   = df['f1_10gb'].rank(method='min', ascending=False).astype(int)  # 1 = most accurate

# y = N - rank so rank 1 sits at the top.
def y_of(rank):
    return N - rank

# Axis definitions: (x position, title, rank series, per-tool value string).
AXES = [
    (0.0, 'Speed',    speed_rank, {t: fmt_time(df.loc[t, 'user_min_10gb']) for t in df.index},
     'single-core CPU time, 10 Gb'),
    (1.0, 'Memory',   mem_rank,   {t: f"{df.loc[t, 'ram_gb_10gb']:.1f} GB" for t in df.index},
     'peak RAM, 10 Gb'),
    (2.0, 'Accuracy', acc_rank,   {t: f"{df.loc[t, 'f1_10gb']:.3f}" for t in df.index},
     'F1, 10 Gb'),
]
XS = [a[0] for a in AXES]

# Per-tool (x, y) across the three axes.
coords = {t: [(x, y_of(r[t])) for (x, _, r, _, _) in AXES] for t in df.index}


# ----- Figure -----
fig, ax = plt.subplots(figsize=(9, 6))
ax.set_xlim(-0.75, 3.02)
ax.set_ylim(-0.85, 6.95)
ax.axis('off')

TOP_Y = N - 1          # y of the best (rank 1) row = 5
BOT_Y = 0              # y of the worst (rank 6) row = 0

# Faint horizontal rank rows spanning the three axes.
for yy in range(N):
    ax.plot([XS[0], XS[-1]], [yy, yy], color=NEUTRAL['grid'], linewidth=0.8,
            zorder=0.3, solid_capstyle='round')

# Subtle spotlight band behind the top (best) row.
ax.add_patch(FancyBboxPatch((XS[0] - 0.16, TOP_Y - 0.34), (XS[-1] - XS[0]) + 0.32, 0.68,
                            boxstyle='round,pad=0.0,rounding_size=0.10',
                            facecolor=NEUTRAL['ink'], alpha=0.045, linewidth=0,
                            zorder=0.4))

# The three vertical parallel axes.
for x in XS:
    ax.plot([x, x], [BOT_Y - 0.42, TOP_Y + 0.42], color=NEUTRAL['rule'],
            linewidth=1.1, zorder=0.6, solid_capstyle='round')

# ----- Tool lines: white casing then coloured line, so crossings read cleanly -----
for tool in TOOL_ORDER:
    xs = [p[0] for p in coords[tool]]
    ys = [p[1] for p in coords[tool]]
    ax.plot(xs, ys, color='white', linewidth=5.4, solid_capstyle='round',
            solid_joinstyle='round', zorder=2)
for tool in TOOL_ORDER:
    xs = [p[0] for p in coords[tool]]
    ys = [p[1] for p in coords[tool]]
    ax.plot(xs, ys, color=TOOL_COLORS[tool], linewidth=2.6, solid_capstyle='round',
            solid_joinstyle='round', alpha=0.92, zorder=3)

# ----- Rank dots -----
for tool in TOOL_ORDER:
    for (x, y) in coords[tool]:
        ax.scatter(x, y, s=76, color=TOOL_COLORS[tool], edgecolor='white',
                   linewidth=1.3, zorder=5)

# ----- End labels: name + value at Speed (left) and Accuracy (right) -----
for tool in df.index:
    c = TOOL_COLORS[tool]
    # Left axis: name + speed value, right-aligned outside the axis.
    (xL, yL) = coords[tool][0]
    ax.annotate(tool, (xL, yL), textcoords='offset points', xytext=(-13, 5.5),
                ha='right', va='center', fontsize=9.2, fontweight='semibold',
                color=c, zorder=6)
    ax.annotate(AXES[0][3][tool], (xL, yL), textcoords='offset points',
                xytext=(-13, -6.5), ha='right', va='center', fontsize=7.6,
                color=NEUTRAL['muted'], zorder=6)
    # Right axis: name + accuracy value, left-aligned outside the axis.
    (xR, yR) = coords[tool][2]
    ax.annotate(tool, (xR, yR), textcoords='offset points', xytext=(13, 5.5),
                ha='left', va='center', fontsize=9.2, fontweight='semibold',
                color=c, zorder=6)
    ax.annotate('F1 ' + AXES[2][3][tool], (xR, yR), textcoords='offset points',
                xytext=(13, -6.5), ha='left', va='center', fontsize=7.6,
                color=NEUTRAL['muted'], zorder=6)

# ----- Middle axis: memory value beside each dot -----
for tool in df.index:
    (xM, yM) = coords[tool][1]
    ax.annotate(AXES[1][3][tool], (xM, yM), textcoords='offset points',
                xytext=(0, 10.5), ha='center', va='bottom', fontsize=7.4,
                color=NEUTRAL['muted'], zorder=6,
                bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                          edgecolor='none', alpha=0.85))

# ----- Axis titles and metric captions -----
for (x, title, _, _, cap) in AXES:
    ax.text(x, 6.62, title, ha='center', va='baseline', fontsize=11.0,
            fontweight='semibold', color=NEUTRAL['ink'])
    ax.text(x, 6.28, cap, ha='center', va='baseline', fontsize=7.9,
            style='italic', color=NEUTRAL['muted'])

# ----- Title -----
ax.set_title('No caller wins on every axis', loc='left', x=-0.11, y=1.0,
             fontsize=12.8, fontweight='semibold', color=NEUTRAL['ink'], pad=14)

fig.tight_layout()
save_figure(fig, 'Fig1B', OUTDIR)
print('done')
