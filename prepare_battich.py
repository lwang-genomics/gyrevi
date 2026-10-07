"""Turn the Battich RPE1-FUCCI scEU-seq data into a GyreVI-ready benchmark.

WHY THIS DATASET. It is the only one here carrying spliced + unspliced *and* an
orthogonal PROTEIN-level phase readout (FUCCI GFP/geminin and RFP/Cdt1 intensities).
Every other reference available to this project -- DeepCycle theta, Seurat scores -- is
itself an RNA-derived estimate, so optimising against them risks reproducing their
assumptions. DeepCycle in particular *assumes* a cycle and assigns every cell to one of
50 bins, so its full circular coverage cannot be used as evidence that our arc is wrong.

WHAT IT SETTLED. The FUCCI angle covers the circle: circ 5-95% = 314 deg, resultant
R = 0.216, 0/12 empty 30-degree sectors, on 5,422 cells. So a cycling population really
does occupy the full circle, and GyreVI's 126 deg / 6-of-12-empty arc on the fibroblasts
is a MODEL FAILURE, not a property of the biology. Note R = 0.216 is *less* concentrated
than DeepCycle's 0.510 -- the truth is closer to uniform than the reference we had been
chasing.

LAYERS. scEU-seq splits reads four ways by splicing and metabolic label:
    sl = spliced labelled     su = spliced unlabelled
    ul = unspliced labelled   uu = unspliced unlabelled
Total spliced = sl + su and total unspliced = ul + uu are independent of labelling
duration, so Pulse and Chase cells at every time point can be pooled -- all 5,422 cells
are usable. (The label split is what scEU-seq is for; GyreVI does not use it.)

Source: Battich et al., Science 2020 (GSE128365), via the figshare copy dynamo ships as
dyn.sample_data.scEU_seq_rpe1().  Download:
    curl -L -o data/battich_rpe1/rpe1.h5ad https://ndownloader.figshare.com/files/47439641
(495 MB. figshare's presigned S3 URL expires in 10 s, so a resumed download appends a
second copy and silently corrupts the file -- always re-download clean, never `curl -C -`.)

Usage:  uv run python prepare_battich.py
"""
import numpy as np, pandas as pd, scanpy as sc, scvelo as scv
import scipy.sparse as sp

from _constants import CC_GENES_HUMAN, get_cc_genes

IN  = "data/battich_rpe1/rpe1.h5ad"
OUT = "data/battich_rpe1/rpe1_gyrevi.h5ad"
MAP = "data/battich_rpe1/ens2sym.tsv"
N_HVG = 2000
# KNN smoothing width. scVelo's default is 30; 15 is measurably better here — see below.
N_NEIGHBORS = 15

a = sc.read_h5ad(IN)
print(f"loaded {a.n_obs} x {a.n_vars}", flush=True)

a.layers["spliced"]   = a.layers["sl"] + a.layers["su"]
a.layers["unspliced"] = a.layers["ul"] + a.layers["uu"]
for k in ("sl", "su", "ul", "uu"):
    del a.layers[k]
a.X = a.layers["spliced"].copy()

# FUCCI reference angle. Each channel is z-scored first: GFP and RFP have different
# dynamic ranges, and atan2 on unscaled channels would bias the angle toward whichever
# channel is larger rather than tracking the cycle.
g = a.obs["GFP_log10_corrected"].values.astype(float)
r = a.obs["RFP_log10_corrected"].values.astype(float)
gz, rz = (g - g.mean()) / g.std(), (r - r.mean()) / r.std()
a.obs["fucci_theta"] = np.arctan2(gz, rz).astype(np.float32)   # (-pi, pi]
a.obs["fucci_magnitude"] = np.sqrt(gz**2 + rz**2).astype(np.float32)
# The dataset ships ENSEMBL ids, so every symbol-keyed gene list (Tirosh CC, the
# selector's readable output) silently matches ZERO genes without this mapping.
# ens2sym.tsv is cached from mygene.info; regenerate it only if the file is missing.
# 11,534 / 11,848 ids map, recovering 91 of the 93 Tirosh genes.
ens2sym = pd.read_csv(MAP, sep="\t", index_col=0)["symbol"].to_dict()
a.var["Gene_Id"] = a.var["Gene_Id"].astype(str)
a.var["symbol"] = [ens2sym.get(g, g) for g in a.var["Gene_Id"]]
a.var_names = a.var["symbol"].values
a.var_names_make_unique()

# Linear-normalised, NOT log1p — GyreVI's decoder reconstructs in linear space and applies
# log1p only inside the loss. Running sc.pp.log1p here would double-log the targets.
# Called as three explicit steps: scvelo 0.3.4's filter_and_normalize forwards
# n_top_genes into normalize_per_cell and raises, and filter_genes_dispersion was
# removed from that release. HVG selection runs BEFORE normalisation because
# seurat_v3 expects raw counts. log1p is deliberately skipped throughout -- GyreVI
# reconstructs in linear space and applies log1p only inside its loss.
scv.pp.filter_genes(a, min_shared_counts=20)
# Never drop a cycling gene: HVG selection on its own discarded every CC gene here,
# which left score_genes_cell_cycle with nothing to score.
sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat_v3", subset=False)
keep = np.asarray(a.var["highly_variable"].values, dtype=bool)
keep |= a.var_names.isin(CC_GENES_HUMAN)
a._inplace_subset_var(keep)
print(f"  kept {a.n_vars} genes ({int(a.var_names.isin(CC_GENES_HUMAN).sum())} CC)", flush=True)
scv.pp.normalize_per_cell(a)
# n_neighbors=15, NOT scVelo's default of 30. Measured on Battich against FUCCI, gene set
# held fixed and the delta target fixed on the unsmoothed layers so neither could drift
# (2 seeds; knn_stage1/2):
#
#   n_neighbors     theta   delta_g    phi_g   discovery      lag fidelity   selector
#   30 (default)   +0.483    +0.309   +0.480       0.931            +0.658      0.854
#   15             +0.501    +0.351   +0.545       0.935            +0.691      0.862
#   10             +0.478    +0.406   +0.514       0.937            +0.722      0.859
#
# 15 beats 30 on every column. delta_g is the one that moves most, and it keeps improving
# down to 10 (+0.406) -- switch to 10 if per-gene kinetics is the headline, at the cost of
# theta and phi_g.
#
# WHY IT HELPS, AND WHY IT IS NOT ABOUT PHASE. A ridge probe of Ms against FUCCI reads
# +0.459 at n=30 and +0.461 at n=10 -- flat -- then jumps to +0.689 with NO smoothing at
# all. Phase loss is a cliff, not a slope: smoothing at all costs a third of the signal
# and the width barely matters (which is why the encoder reads `rawln` instead). What the
# width does affect is the u/s LAG: agreement with the unsmoothed lag runs +0.658 at n=30
# to +0.765 at n=5, and delta_g is fitted against exactly that.
#
# Do NOT drop smoothing entirely. The selector needs it: |corr(s,u)| on cycling genes
# falls 0.481 -> 0.085 unsmoothed and its AUROC drops 0.854 -> 0.780. The coherence the
# selector keys on is largely manufactured by smoothing, amplifying a weak real signal.
scv.pp.moments(a, n_pcs=30, n_neighbors=N_NEIGHBORS)
# Seurat S/G2M scores. GyreVI.setup_anndata() looks for these and silently disables
# L_align if they are absent (printing a warning that is easy to miss in a long log),
# so they are computed here rather than left to the caller. Scored on log1p of a COPY:
# score_genes_cell_cycle expects log data, while the saved layers must stay linear.
cc_genes, _ = get_cc_genes("human")
s_genes   = [g for g in CC_GENES_HUMAN[:43] if g in a.var_names]     # G1/S block
g2m_genes = [g for g in CC_GENES_HUMAN[43:] if g in a.var_names]     # G2/M block
tmp = a.copy(); sc.pp.log1p(tmp)
sc.tl.score_genes_cell_cycle(tmp, s_genes=s_genes, g2m_genes=g2m_genes)
for k in ("S_score", "G2M_score", "phase"):
    a.obs[k] = tmp.obs[k].values
print(f"  scored cell cycle on {len(s_genes)} S + {len(g2m_genes)} G2M genes present; "
      f"phases {dict(a.obs['phase'].astype(str).value_counts())}", flush=True)

# Unsmoothed encoder input. scv.pp.moments averages each cell with its neighbours at
# similar-but-different phases, which is a moving average over a periodic signal: an
# identical ridge probe reaches +0.462 against FUCCI on Ms but +0.688 on the same cells
# unsmoothed. The decoder still needs the smoothed moments for velocity kinetics, so this
# layer feeds the theta encoder only, via GyreVI.setup_anndata(theta_layer="rawln"),
# and the z encoder via z_layer="rawln" (both are setup_anndata defaults).
_raw = np.asarray(a.layers["spliced"].todense() if sp.issparse(a.layers["spliced"])
                  else a.layers["spliced"], dtype=np.float32)
_sz = _raw.sum(1, keepdims=True); _sz[_sz == 0] = 1.0
a.layers["rawln"] = np.log1p(_raw / _sz * np.median(_raw.sum(1))).astype(np.float32)

# Same transform applied to UNSPLICED. theta_layer only ever replaced the SPLICED channel,
# so before this the theta encoder read unsmoothed log spliced against smoothed LINEAR Mu --
# mismatched in both scale and smoothing. Measured on Battich, 4 seeds, paired within one
# process, against the previous default (theta rawln + Mu, z Ms):
#
#   theta encoder unspliced channel      theta    delta_g     phi_g   discovery
#   Mu            (linear, smoothed)    +0.473     +0.592    +0.684       0.909
#   log1p(Mu)     (log,    smoothed)    +0.461     +0.592    +0.681       0.942
#   rawln_u       (log,  unsmoothed)    +0.467     +0.596    +0.690       0.939
#
# The LOG is what moves discovery: both logged arms gain ~+0.030 (4/4 seeds) while
# discovery is otherwise pinned at 0.908-0.910 across every fit with seed sd ~0.001.
# The UNSMOOTHING then recovers the kinetics: delta_g and phi_g go from -0.002/-0.003
# under log1p(Mu) to +0.005/+0.006 (both 4/4) under rawln_u. Cost is theta -0.006 against
# FUCCI (3/4 worse) -- note phi_g uses a rotation-invariant circular correlation, so the
# gene-phase PATTERN improving while absolute FUCCI agreement dips is self-consistent.
# 86.6% of raw unspliced entries are zero (vs 64.1% spliced); it still wins.
_rawu = np.asarray(a.layers["unspliced"].todense() if sp.issparse(a.layers["unspliced"])
                   else a.layers["unspliced"], dtype=np.float32)
_szu = _rawu.sum(1, keepdims=True); _szu[_szu == 0] = 1.0
a.layers["rawln_u"] = np.log1p(_rawu / _szu * np.median(_rawu.sum(1))).astype(np.float32)

a.write(OUT)
print(f"wrote {OUT}  {a.n_obs} x {a.n_vars}  layers {sorted(a.layers)}", flush=True)

t = a.obs["fucci_theta"].values
mu = np.angle(np.mean(np.exp(1j*t))); tc = np.angle(np.exp(1j*(t-mu)))
q = np.degrees(np.quantile(tc, [.05, .95]))
h = np.histogram(tc, bins=12, range=(-np.pi, np.pi))[0]
print(f"FUCCI reference: circ5-95 {q[1]-q[0]:.0f}deg  R {abs(np.mean(np.exp(1j*t))):.3f}  "
      f"empty {int((h==0).sum())}/12   <- the target GyreVI must reproduce")
