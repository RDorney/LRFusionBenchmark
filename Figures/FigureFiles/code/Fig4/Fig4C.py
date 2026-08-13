"""Figure 4C. Known-fusion recovery in Huh7 across libraries.

Each column is one literature-known fusion, grouped into HCC and Huh7 source
categories (colour stripe at the bottom). Each row is a library type. A filled
dot means at least one tool or replicate of that library detected the fusion, by
Ensembl gene-ID match; dots in a column are joined to show co-detection. Only
fusions detected in at least one library are shown. The right-hand bars give the
number of known fusions each library recovered.

Data preparation and gene-ID matching live in _known_fusion_data.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
from _paths import OUT
OUTDIR = OUT / 'Fig4'
from palette import NEUTRAL, set_style, save_figure
from _known_fusion_data import prepare, LIBRARY_ROWS

set_style()

# reads >= 2 or pairs >= 2: the analysis-wide definition of a fusion call, matching
# every other display item. A no-op here (curated known fusions all clear it), applied
# so the read-support filter is uniform across the paper.
all_fusions, category_idx, M = prepare(read_support=2, restrict_known=False,
                                       detected_only=True)
n_cols = len(all_fusions)
n_rows = len(LIBRARY_ROWS)
n_det = M.sum(axis=0)
lib_recov = M.sum(axis=1)

CAT_COLORS = {'HCC': '#3973A0', 'HUH7': '#5AAE61'}
CAT_LABEL = {'HCC': 'HCC', 'HUH7': 'Huh7'}


# ---------- Figure ----------
fig = plt.figure(figsize=(12.0, 6.8))
fig.patch.set_facecolor('white')

ML, MR = 0.150, 0.850          # matrix left / right (fraction)
MW = MR - ML

CAT_B, CAT_H = 0.880, 0.030    # category stripe (top)
DOT_B, DOT_H = 0.610, 0.250    # dot grid
NAMES_B, NAMES_H = 0.095, 0.495  # fusion names (bottom)
BAR_L = MR + 0.020
BAR_W = 0.965 - BAR_L

# column guides from the category stripe down to the dot grid
axg = fig.add_axes([ML, NAMES_B + NAMES_H, MW, CAT_B - (NAMES_B + NAMES_H)])
axg.set_xlim(-0.5, n_cols - 0.5)
axg.set_ylim(0, 1)
axg.set_axis_off()
for j in range(n_cols):
    axg.plot([j, j], [0, 1], color='#ECEFF1', linewidth=0.5, zorder=0)


# ----- dot grid -----
ax = fig.add_axes([ML, DOT_B, MW, DOT_H])
ax.set_xlim(-0.5, n_cols - 0.5)
ax.set_ylim(-0.5, n_rows - 0.5)
ax.invert_yaxis()
ax.set_axis_off()

for j in range(n_cols):
    color = CAT_COLORS[category_idx[j]]
    present = np.where(M[:, j])[0]
    if len(present) >= 2:
        ax.plot([j, j], [present.min(), present.max()], color=color,
                linewidth=1.6, zorder=3, solid_capstyle='round')
    for i in range(n_rows):
        if M[i, j]:
            ax.scatter(j, i, s=95, c=color, edgecolors='white',
                       linewidths=0.5, zorder=4)
        else:
            ax.scatter(j, i, s=9, c='#D6DBE0', edgecolors='none', zorder=2)

# library row labels (left of grid)
for i, (_, _, label) in enumerate(LIBRARY_ROWS):
    ax.text(-0.9, i, label, ha='right', va='center', fontsize=9.4,
            color=NEUTRAL['ink'], fontweight='semibold', clip_on=False)


# ----- right-hand recovery bars -----
axb = fig.add_axes([BAR_L, DOT_B, BAR_W, DOT_H])
axb.set_xlim(0, n_cols)
axb.set_ylim(-0.5, n_rows - 0.5)
axb.invert_yaxis()
axb.set_axis_off()
for i in range(n_rows):
    axb.barh(i, lib_recov[i], height=0.46, color=NEUTRAL['muted'],
             edgecolor='none', zorder=2)
    axb.text(lib_recov[i] + n_cols * 0.02, i, f'{lib_recov[i]}',
             ha='left', va='center', fontsize=8.8, fontweight='semibold',
             color=NEUTRAL['ink'])
axb.text(0, -0.85, f'Known recovered  (of {n_cols})', ha='left', va='bottom',
         fontsize=8.2, color=NEUTRAL['muted'])


# ----- fusion name labels (below grid, reading down) -----
axn = fig.add_axes([ML, NAMES_B, MW, NAMES_H])
axn.set_xlim(-0.5, n_cols - 0.5)
axn.set_ylim(0, 1)
axn.set_axis_off()
for j, fname in enumerate(all_fusions):
    axn.text(j, 0.99, fname, ha='center', va='top', rotation=90,
             fontsize=6.7, color=NEUTRAL['body'])


# ----- category stripe (bottom) -----
axc = fig.add_axes([ML, CAT_B, MW, CAT_H])
axc.set_xlim(-0.5, n_cols - 0.5)
axc.set_ylim(0, 1)
axc.set_axis_off()
start = 0
for cat in ['HCC', 'HUH7']:
    n = sum(1 for c in category_idx if c == cat)
    if n == 0:
        continue
    axc.add_patch(Rectangle((start - 0.5, 0), n, 1, facecolor=CAT_COLORS[cat],
                            edgecolor='white', linewidth=1.2, zorder=2))
    axc.text(start + n / 2 - 0.5, 0.5, CAT_LABEL[cat], ha='center', va='center',
             fontsize=9.6, fontweight='bold', color='white', zorder=3)
    start += n


# ---------- Title ----------
fig.text(0.035, 0.955, 'Known fusion recovery in Huh7 across libraries',
         ha='left', va='top', fontsize=12.8, fontweight='semibold',
         color=NEUTRAL['ink'])


save_figure(fig, 'Fig4C', OUTDIR)
plt.close(fig)
