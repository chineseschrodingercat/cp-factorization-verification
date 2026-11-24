"""Draw Figure 1: the smallest collapse example (g = q = 3, Appendix A).

(a) Before deletion every cross-group pair of active labels is an edge of the
complete tripartite graph and needs its own column (27 columns, plus the
empty word).  (b) After deleting the slack coordinate z, the nine words
(a, b, a+b) over F_3 form an orthogonal array: as a Latin square, each cell,
each (row, symbol) and each (column, symbol) pair occurs exactly once.
Run:  python make_figure_collapse.py   (writes fig_collapse.pdf)
"""
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})

fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.9), gridspec_kw={"width_ratios": [1.15, 1]})

# (a) the complete tripartite graph of cross-group pairs
ax = axes[0]
centers = {0: np.array([0.0, 0.0]), 1: np.array([4.0, 0.0]), 2: np.array([2.0, 3.4])}
centroid = sum(centers.values()) / 3
pos = {}
for g, c in centers.items():
    d = centroid - c
    u = np.array([-d[1], d[0]]) / np.linalg.norm(d)
    for a in range(3):
        pos[(g, a)] = c + (a - 1) * 0.55 * u
for (g1, a), (g2, b) in itertools.combinations(pos, 2):
    if g1 != g2:
        ax.plot(*zip(pos[(g1, a)], pos[(g2, b)]), color="0.6", lw=0.9, zorder=1)
for (g, a), p in pos.items():
    ax.plot(*p, "o", ms=5.5, color="black", zorder=5)
    off = (p - centroid) / np.linalg.norm(p - centroid) * 0.36
    ax.text(*(p + off), str(a), ha="center", va="center", fontsize=8.5)
for g, c in centers.items():
    off = (c - centroid) / np.linalg.norm(c - centroid) * 1.05
    ax.text(*(c + off), rf"$G_{g + 1}$", ha="center", va="center", fontsize=11)
ax.set_aspect("equal"); ax.axis("off")
ax.set_xlim(-1.3, 5.3); ax.set_ylim(-1.15, 4.75)
ax.set_title("(a) with $z$: one column per edge", fontsize=9.5)

# (b) the orthogonal array (a, b, a+b) as a Latin square
ax = axes[1]
colors = [plt.cm.tab10(i) for i in (0, 1, 2, 3, 4, 5, 6, 8, 9)]
for k, (a, b) in enumerate(itertools.product(range(3), repeat=2)):
    x, y = b, 2 - a
    ax.add_patch(Rectangle((x, y), 1, 1, facecolor=colors[k], alpha=0.30, edgecolor="black", lw=0.8))
    ax.text(x + 0.5, y + 0.5, str((a + b) % 3), ha="center", va="center", fontsize=13)
for a in range(3):
    ax.text(-0.25, 2.5 - a, str(a), ha="center", va="center", fontsize=9)
    ax.text(a + 0.5, 3.25, str(a), ha="center", va="center", fontsize=9)
ax.text(-0.75, 1.5, "$G_1$", ha="center", va="center", fontsize=11)
ax.text(1.5, 3.7, "$G_2$", ha="center", va="center", fontsize=11)
ax.text(3.45, 1.5, "entry\nin $G_3$", ha="center", va="center", fontsize=9)
ax.set_aspect("equal"); ax.axis("off")
ax.set_xlim(-1.1, 4.1); ax.set_ylim(-0.35, 4.3)
ax.set_title("(b) without $z$: nine words $(a,b,a{+}b)$", fontsize=9.5)

fig.tight_layout(w_pad=1.0)
fig.savefig("fig_collapse.pdf", bbox_inches="tight")
print("wrote fig_collapse.pdf")
