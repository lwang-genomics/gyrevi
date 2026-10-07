"""scVI baseline: how much cell cycle phase does a standard latent carry?

Trains vanilla scVI on the same cells/genes and linearly probes its latent for phase.
This is the "what does a generic deep generative model give you for free" bar, alongside
Seurat (prior knowledge, no fitting) and DeepCycle (dedicated cycle method).

The probe is SUPERVISED (5-fold CV ridge to cos/sin of the reference), so it measures
information PRESENT in the latent, not what scVI would report unsupervised. That makes it
an upper bound on scVI's phase content and a fair ceiling to compare theta against.
"""
import sys, numpy as np, scanpy as sc, scvi, torch
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import cross_val_predict
from circular import compare_phases, seurat_angle, to_radians, wrap_2pi

scvi.settings.seed = 0
DS = sys.argv[1] if len(sys.argv) > 1 else "fibroblast"

if DS == "battich":
    a = sc.read_h5ad("data/battich_rpe1/rpe1_gyrevi.h5ad")
    counts_layer = "spliced"
    refs = {"FUCCI": wrap_2pi(a.obs["fucci_theta"].values.astype(float))}
else:
    a = sc.read_h5ad("data/fibroblast/velocity_anndata_human_fibroblast_DeepCycle_ISMARA.h5ad")
    d = np.load("data/fibroblast/theta_diag_2004.npz", allow_pickle=True)
    a = a[:, np.isin(a.var_names, np.asarray(d["var_names"]))].copy()
    counts_layer = "spliced"
    refs = {"DeepCycle": wrap_2pi(to_radians(a.obs["cell_cycle_theta"].values.astype(float), 1.0))}
refs["Seurat"] = wrap_2pi(seurat_angle(a.obs["S_score"].values, a.obs["G2M_score"].values))
print(f"{DS}: {a.n_obs} cells x {a.n_vars} genes", flush=True)

X = a.layers[counts_layer]
a.layers["counts"] = np.rint(np.asarray(X.todense() if hasattr(X, "todense") else X)).astype(np.float32)
scvi.model.SCVI.setup_anndata(a, layer="counts")
model = scvi.model.SCVI(a, n_latent=10, n_hidden=128, n_layers=1)
model.train(max_epochs=300, early_stopping=False, accelerator="mps",
            plan_kwargs={"lr": 1e-3})
Z = model.get_latent_representation()
print(f"scVI latent {Z.shape}", flush=True)

print(f"\n{'reference':<14}{'scVI probe':>12}{'Seurat':>10}{'medErr':>9}")
print("-"*46)
for name, r in refs.items():
    tgt = np.c_[np.cos(r), np.sin(r)]
    p = cross_val_predict(RidgeCV(alphas=np.logspace(-2, 4, 14)), Z, tgt, cv=5)
    c = compare_phases(r, wrap_2pi(np.arctan2(p[:, 1], p[:, 0])), name, "scVI")
    cs = compare_phases(r, refs["Seurat"], name, "Seurat")
    print(f"{name:<14}{c['spearman']:>+12.3f}{cs['spearman']:>+10.3f}"
          f"{c['median_abs_err_deg']:>8.0f}°")
import os
np.save(os.path.join(os.environ.get("SCR_OUT", "."), f"scvi_latent_{DS}.npy"), Z)
print("saved latent")
