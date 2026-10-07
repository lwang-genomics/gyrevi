"""Run CycleVI (Mozdzanowski et al. 2025) on Battich RPE1-FUCCI.

Their paper validates on THIS dataset and reports Spearman 0.513 against FUCCI (Fig. 5),
so this is a direct head-to-head. Running it ourselves rather than quoting that number
means identical preprocessing, identical cells and our own circular metric.

CycleVI is spliced-only by design; phase is the angular coordinate of the first two
latent dimensions. Initialisation follows their sec 2.1.3: the Seurat atan2(G2M, S) angle,
quantile-transformed to uniform on [0, 2pi).
"""
import numpy as np, scanpy as sc, scipy.sparse as sp, pandas as pd, torch
from pathlib import Path
from cyclevi import CycleVI, create_cell_cycle_gene_mask, get_cc_genes_path

OUT = str(Path(__file__).resolve().parent.parent / "data/battich_rpe1/cyclevi_theta.csv")
prep = sc.read_h5ad("data/battich_rpe1/rpe1_gyrevi.h5ad")   # for gene subset + obs
raw  = sc.read_h5ad("data/battich_rpe1/rpe1.h5ad")

# NB likelihood needs counts: rebuild total spliced from the scEU-seq split, before the
# per-cell normalisation prepare_battich.py applies.
raw.layers["spliced_counts"] = raw.layers["sl"] + raw.layers["su"]
raw.var["Gene_Id"] = raw.var["Gene_Id"].astype(str)
keep = raw.var["Gene_Id"].isin(set(prep.var["Gene_Id"].astype(str)))
a = raw[:, keep.values].copy()
assert list(a.obs_names) == list(prep.obs_names), "cell order mismatch"
X = a.layers["spliced_counts"]
a.layers["counts"] = np.rint(np.asarray(X.todense() if sp.issparse(X) else X)).astype(np.float32)
for k in ("S_score", "G2M_score", "fucci_theta"):
    a.obs[k] = prep.obs[k].values
print(f"{a.n_obs} cells x {a.n_vars} genes (counts)", flush=True)

# Seurat init angle, quantile-transformed to uniform (their sec 2.1.3).
raw_ang = np.arctan2(a.obs["G2M_score"].values.astype(float),
                     a.obs["S_score"].values.astype(float))
o = raw_ang.argsort(kind="mergesort"); rk = np.empty(a.n_obs); rk[o] = np.arange(a.n_obs)
a.obs["cycle_init_angle"] = ((rk / a.n_obs) * 2 * np.pi).astype(np.float32)

# NOTE: cyclevi's own create_cell_cycle_gene_mask() cannot read the file cyclevi ships.
# It does pd.read_csv(path, header=None)[0], i.e. column 0, but the bundled
# homo_sapiens_cc_genes.csv has columns phase,geneID,modified -- so column 0 is the PHASE
# label ("G2/M"), not the gene id, and the mask comes back with ZERO genes. Running CycleVI
# that way silently disables its cell-cycle gene mask and handicaps it, so the mask is
# built here from the geneID column instead. Verified non-empty below.
cc = pd.read_csv(get_cc_genes_path())
ids = set(cc["geneID"].astype(str).str.upper())
mask = torch.tensor([g.upper() in ids for g in a.var["Gene_Id"].astype(str)], dtype=torch.bool)
print(f"cycle gene mask: {int(mask.sum())} genes  (their helper returns 0 -- see comment)",
      flush=True)
assert int(mask.sum()) > 50, "mask still empty -- do not report a handicapped baseline"

CycleVI.setup_anndata(a, layer="counts", cycle_initiation_angle_key="cycle_init_angle")
model = CycleVI(a, n_latent=10, n_hidden=128, n_layers=1, cycle_gene_mask=mask)
model.train(max_epochs=400, batch_size=128, early_stopping=False)

z = model.get_latent_representation()
theta = np.arctan2(z[:, 1], z[:, 0]) % (2 * np.pi)
pd.DataFrame({"cyclevi_theta": theta}, index=a.obs_names).to_csv(OUT)
print(f"wrote {OUT}  latent {z.shape}", flush=True)
