#!/usr/bin/env python3
"""Draw Figure 2 of the paper from benchmark/results/summary.tsv (requires matplotlib)."""

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results", "summary.tsv")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results", "figure2")

rows = list(csv.DictReader(open(src), delimiter="\t"))
scen = list(dict.fromkeys(r["scenario"] for r in rows))
get = {(r["scenario"], r["condition"]): r for r in rows}
style = {
    "error-free": dict(color="#9a9a93", marker="o", mfc="#9a9a93", off=-0.22, label="Error-free"),
    "noisy": dict(color="#2a6fdb", marker="o", mfc="#2a6fdb", off=0.0, label="Noisy, strict"),
    "noisy, missing ignored": dict(color="#2a6fdb", marker="o", mfc="white", off=0.22,
                                   label="Noisy, missing genotypes ignored"),
}

plt.rcParams.update({"font.size": 8, "font.family": "sans-serif", "axes.spines.top": False,
                     "axes.spines.right": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.3), sharey=True,
                           gridspec_kw={"width_ratios": [1.5, 1]})
ys = range(len(scen))
for cond, st in style.items():
    a.plot([float(get[(s, cond)]["sensitivity"]) * 100 for s in scen],
           [y + st["off"] for y in ys], ls="none", marker=st["marker"], ms=5,
           mec=st["color"], mfc=st["mfc"], mew=1.3, label=st["label"])
a.set_xlim(80, 101)
a.set_xticks([80, 85, 90, 95, 100])
a.set_xlabel("Families with causal variant retained (%)")
a.set_yticks(list(ys))
a.set_yticklabels(scen)
a.invert_yaxis()
a.grid(axis="x", color="#e6e6e6", lw=0.6)
a.set_axisbelow(True)
fig.legend(*a.get_legend_handles_labels(), loc="upper center", ncol=3, frameon=False,
           fontsize=7, bbox_to_anchor=(0.5, 1.0))
a.set_title("A", loc="left", fontweight="bold")

for cond in ("error-free", "noisy"):
    st = style[cond]
    vals = [float(get[(s, cond)]["median_background"]) for s in scen]
    yy = [y + (-0.18 if cond == "error-free" else 0.18) for y in ys]
    b.barh(yy, vals, height=0.34, color=st["color"])
    for v, y in zip(vals, yy):
        b.text(v + 4, y, ("%g" % v), va="center", fontsize=6.5, color="#444")
b.set_xlabel("Median other variants reported")
b.set_xlim(0, 340)
b.set_title("B", loc="left", fontweight="bold")

fig.tight_layout(rect=(0, 0, 1, 0.93))
for ext in ("pdf", "png"):
    fig.savefig("%s.%s" % (out, ext), dpi=300)
print("written %s.pdf and %s.png" % (out, out))
