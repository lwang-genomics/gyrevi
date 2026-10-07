"""Phase-recovery comparison against the FUCCI protein ground truth.

Self-contained: colours are inlined rather than imported from a local style file, so
this reproduces anywhere. Numbers are the 6-seed paired benchmark.

Two conventions matter and the figure states both:
  * the oracle row is FITTED to the ground-truth labels, so it is a ceiling on the phase
    information present in the data, not a competing method;
  * CycleVI's score depends on the alignment convention (circular alignment vs a
    best-shift search), so the bar carries its range rather than one chosen number.
"""
import matplotlib.pyplot as plt

SIG, NONSIG, RULE, LIGHT = "#D95F02", "#BDBDBD", "#2166AC", "#E8E8E8"
plt.rcParams.update({"font.size": 14, "figure.facecolor": "white", "axes.facecolor": "white"})

# label, score, upper range (None = single-valued), is_oracle
ROWS = [
    ("DeepCycle",                     0.275, None,  False),
    ("CycleVI",                       0.438, 0.573, False),
    ("Seurat angle\n(two gene-set means)", 0.463, None,  False),
    ("GyreVI $\\theta$",              0.473, None,  False),
    ("scVI latent + ridge probe",     0.669, None,  True),
]

fig, ax = plt.subplots(figsize=(10, 5))
y = range(len(ROWS))

for i, (label, val, hi, oracle) in enumerate(ROWS):
    if hi is not None:                       # convention-dependent range
        ax.barh(i, hi - val, left=val, color=LIGHT, edgecolor="#999999",
                linewidth=0.7, height=0.65, hatch="///", zorder=2)
    colour = SIG if label.startswith("GyreVI") else NONSIG
    ax.barh(i, val, color=colour, edgecolor="#4D4D4D", linewidth=0.7, height=0.65, zorder=3)
    if oracle:
        ax.text(val - 0.012, i, "oracle — fitted to labels", ha="right", va="center",
                fontsize=12, color="#4D4D4D", zorder=4)

ax.axvline(0.463, color=RULE, linestyle="--", linewidth=1.5, zorder=1)
ax.text(0.463, len(ROWS) - 0.35, " label-free baseline", color=RULE, fontsize=13,
        ha="left", va="bottom")
ax.text(0.585, 1.0, "range spans the\nalignment convention", fontsize=11,
        color="#6E6E6E", ha="left", va="center")

ax.set_yticks(list(y))
ax.set_yticklabels([r[0] for r in ROWS])
ax.set_xlabel("Spearman correlation with FUCCI protein phase")
ax.set_xlim(0, 0.80)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.grid(False)
plt.tight_layout()
fig.savefig("phase_method_comparison.png", dpi=300, bbox_inches="tight")
print("wrote figures/phase_method_comparison.png")
