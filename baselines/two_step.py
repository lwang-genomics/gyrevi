"""The simple estimators that beat GyreVI — reproducible from this repository alone.

A phase from a label-free method, then the same least-squares cosine fit per gene that
defines the per-gene targets. Scored exactly as the GyreVI benchmark scores the model:
the targets are cosine fits at the PROTEIN phase, on the same smoothed counts.

    python baselines/two_step.py [path/to/rpe1_gyrevi.h5ad]     # made by prepare_battich.py
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from _constants import CC_GENES_HUMAN
from circular import align_circular, circ_corr, compare_phases, seurat_angle, wrap_2pi, wrap_pi
from gene_selection import cycling_scores

path = sys.argv[1] if len(sys.argv) > 1 else ROOT / "data/battich_rpe1/rpe1_gyrevi.h5ad"
a = ad.read_h5ad(path)
dense = lambda X: np.asarray(X.todense() if sp.issparse(X) else X, dtype=np.float64)
Ms, Mu = dense(a.layers["Ms"]), dense(a.layers["Mu"])
protein = wrap_2pi(a.obs["fucci_theta"].values.astype(float))

# ── two label-free phases ─────────────────────────────────────────────────────────
seurat = wrap_2pi(seurat_angle(a.obs["S_score"].values, a.obs["G2M_score"].values))
cc = np.asarray([g in set(CC_GENES_HUMAN) for g in a.var_names])
Y = np.log1p(dense(a.layers["spliced"])[:, cc])
Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-8)
P = PCA(n_components=2, whiten=True, random_state=0).fit_transform(Y)
pca = wrap_2pi(np.arctan2(P[:, 1], P[:, 0]))

# ── per-gene targets: a cosine fit at the protein phase ──────────────────────────
def cosine_fit(phase):
    D = np.c_[np.ones_like(phase), np.cos(phase), np.sin(phase)]
    cs, *_ = np.linalg.lstsq(D, np.log1p(Ms), rcond=None)
    cu, *_ = np.linalg.lstsq(D, np.log1p(Mu), rcond=None)
    cl, *_ = np.linalg.lstsq(D, Ms, rcond=None)
    return dict(peak_s=np.arctan2(cs[2], cs[1]), peak_u=np.arctan2(cu[2], cu[1]),
                amp_s=np.hypot(cs[1], cs[2]), amp_u=np.hypot(cu[1], cu[2]), amp=np.hypot(cl[1], cl[2]))

T = cosine_fit(protein)
osc = np.zeros(a.n_vars, bool)
osc[np.argsort(cycling_scores(Ms, Mu)["combined"])[::-1][:100]] = True      # model-independent
strong = osc & (T["amp_s"] > np.quantile(T["amp_s"][osc], 0.25))
fair = strong.copy()                    # the lead is only estimable where unspliced oscillates
fair[strong] = T["amp_u"][strong] >= np.median(T["amp_u"][strong])
lead = wrap_pi(T["peak_s"] - T["peak_u"])

print(f"{a.n_obs} cells | {cc.sum()} curated cell-cycle genes | "
      f"{osc.sum()} oscillating / {strong.sum()} strong / {fair.sum()} unspliced-informative genes\n")
print(f"{'phase from':12s} {'cell phase':>11s} {'peak phase':>11s} {'lead':>8s} {'amplitude':>10s}")
for name, phase in [("Seurat", seurat), ("PCA", pca)]:
    theta = compare_phases(protein, phase)["spearman"]
    f = cosine_fit(align_circular(protein, phase)[0])     # direction + rotation only
    phi = circ_corr(wrap_2pi(T["peak_s"][strong]), wrap_2pi(f["peak_s"][strong]))
    delta = spearmanr(wrap_pi(f["peak_s"] - f["peak_u"])[fair], lead[fair]).statistic
    R = spearmanr(f["amp"][osc], T["amp"][osc]).statistic
    print(f"{name:12s} {theta:>+11.3f} {phi:>+11.3f} {delta:>+8.3f} {R:>+10.3f}")
print(f"{'GyreVI':12s} {'+0.473':>11s} {'+0.706':>11s} {'+0.598':>8s} {'+0.912':>10s}   (learned, 6-seed mean)")
