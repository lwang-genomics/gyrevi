"""Preprocessing utilities for GyreVI.

GyreVI takes KNN-smoothed moments (Ms, Mu) as input — the same preprocessing
pipeline as veloVI and scVelo dynamical model. This module provides a single
convenience function that wraps the required scVelo + Scanpy steps.

Typical call
------------
>>> import scvelo as scv
>>> import scanpy as sc
>>> from preprocess import preprocess_for_gyrevi
>>>
>>> adata = scv.datasets.pancreas()
>>> preprocess_for_gyrevi(adata, n_top_genes=2000, organism="mouse")
>>> # adata now has layers: Ms, Mu  — ready for GyreVI.setup_anndata()
"""

from __future__ import annotations

from typing import Optional

import scanpy as sc
import scvelo as scv
from anndata import AnnData


def preprocess_for_gyrevi(
    adata: AnnData,
    n_top_genes:    int   = 2000,
    n_pcs:          int   = 30,
    n_neighbors:    int   = 30,
    min_shared_counts: int = 20,
    organism:       str   = "mouse",
    copy:           bool  = False,
    verbose:        bool  = True,
) -> AnnData:
    """Preprocess AnnData for GyreVI input.

    Runs the standard scVelo moment-computation pipeline:

    1. Filter genes with low shared spliced/unspliced counts.
    2. Normalise spliced, unspliced, and X per cell.
    3. Log1p transform X.
    4. Select top highly variable genes (HVGs).
    5. Compute PCA + KNN graph.
    6. Compute KNN-smoothed moments Ms and Mu.

    After this function, ``adata.layers`` contains ``'Ms'`` and ``'Mu'``,
    which are the inputs expected by :meth:`~GyreVI.setup_anndata`.

    Parameters
    ----------
    adata : AnnData
        Raw AnnData with ``layers['spliced']`` and ``layers['unspliced']``.
        Typically loaded via ``scv.datasets.*`` or from a loom/h5ad file with
        scVelo-compatible layers.
    n_top_genes : int
        Number of highly variable genes to retain (default 2000).
        Use more genes for better cycling gene discovery, fewer for speed.
    n_pcs : int
        Number of PCs for KNN graph construction (default 30).
    n_neighbors : int
        Number of neighbours for KNN graph and moment smoothing (default 30).
    min_shared_counts : int
        Minimum total counts in spliced + unspliced for a gene to be retained
        (default 20). Lower = more genes retained; raise if data is sparse.
    organism : str
        ``'human'`` or ``'mouse'`` — used only for informational printing.
    copy : bool
        If ``True``, return a copy and leave the original unchanged.
    verbose : bool
        Print progress messages.

    Returns
    -------
    AnnData
        Preprocessed AnnData with Ms and Mu layers, PCA, UMAP, and neighbours.
        Returns the original object in-place (or a copy if ``copy=True``).

    Notes
    -----
    - **Do not** apply ``sc.pp.log1p`` to Ms/Mu — these are moment-smoothed
      normalised counts and must remain in linear space for the cosinor decoder.
    - The log1p applied here is only to ``adata.X`` for HVG selection; it does
      not touch the Ms/Mu layers.
    - If your data already has Ms and Mu layers (e.g. from a pre-processed
      h5ad file), you can skip this function entirely.

    Examples
    --------
    >>> adata = scv.datasets.pancreas()
    >>> preprocess_for_gyrevi(adata, n_top_genes=2000, organism="mouse")
    >>> print(list(adata.layers.keys()))   # ['spliced', 'unspliced', 'Ms', 'Mu']
    """
    if copy:
        adata = adata.copy()

    if verbose:
        print(f"[preprocess] Input: {adata.n_obs} cells × {adata.n_vars} genes  ({organism})")

    # ------------------------------------------------------------------
    # Step 1 — Filter low-count genes + normalise per cell
    # scv.pp.filter_and_normalize handles both spliced and unspliced
    # ------------------------------------------------------------------
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=min_shared_counts,
        # log=False not passed — not supported in scvelo 0.3.x
        # log1p is applied manually below, after HVG selection
    )
    if verbose:
        print(f"[preprocess] After gene filter: {adata.n_vars} genes")

    # ------------------------------------------------------------------
    # Step 2 — Log1p X for HVG selection (does NOT affect Ms/Mu)
    # ------------------------------------------------------------------
    sc.pp.log1p(adata)

    # ------------------------------------------------------------------
    # Step 3 — Highly variable genes
    # seurat_v3 (VST) is preferred — requires scikit-misc.
    # Falls back to seurat flavor (no extra deps) if skmisc is missing.
    # ------------------------------------------------------------------
    try:
        from skmisc.loess import loess  # noqa: F401
        hvg_flavor = "seurat_v3"
        hvg_layer  = "spliced"          # VST works on normalised counts
    except ImportError:
        if verbose:
            print("[preprocess] scikit-misc not found — falling back to 'seurat' HVG flavor. "
                  "Run `uv add scikit-misc` for the preferred seurat_v3 method.")
        hvg_flavor = "seurat"
        hvg_layer  = None               # seurat flavor uses log-normalised adata.X

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        subset=True,
        layer=hvg_layer,
        flavor=hvg_flavor,
    )
    if verbose:
        print(f"[preprocess] After HVG selection: {adata.n_vars} genes")

    # ------------------------------------------------------------------
    # Step 4 — PCA + KNN neighbours
    # scv.pp.moments also runs PCA/KNN internally, but we do it explicitly
    # so the parameters are visible and stored properly in adata.
    # ------------------------------------------------------------------
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # ------------------------------------------------------------------
    # Step 5 — KNN-smoothed moments: Ms and Mu
    # This is the key step — GyreVI operates on Ms and Mu, not raw counts.
    # ------------------------------------------------------------------
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    if verbose:
        print(f"[preprocess] Computed moments: Ms {adata.layers['Ms'].shape}, "
              f"Mu {adata.layers['Mu'].shape}")

    # ------------------------------------------------------------------
    # Step 6 — UMAP (computed now so it's ready for visualisation)
    # ------------------------------------------------------------------
    sc.tl.umap(adata)
    if verbose:
        print(f"[preprocess] UMAP computed and stored in adata.obsm['X_umap']")
        print(f"[preprocess] Done. Layers available: {list(adata.layers.keys())}")
        print(f"[preprocess] → Call GyreVI.setup_anndata(adata, "
              f"spliced_layer='Ms', unspliced_layer='Mu')")

    return adata


def check_preprocessing(adata: AnnData) -> bool:
    """Verify that adata is ready for GyreVI.

    Checks that required layers and basic data quality conditions are met.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData to check.

    Returns
    -------
    bool
        ``True`` if all checks pass, ``False`` otherwise (with printed warnings).
    """
    ok = True

    # Check layers
    for layer in ("Ms", "Mu"):
        if layer not in adata.layers:
            print(f"[check] ✗ Missing layer '{layer}' — run preprocess_for_gyrevi() first.")
            ok = False
        else:
            import numpy as np
            from scipy import sparse
            data = adata.layers[layer]
            if sparse.issparse(data):
                data = data.toarray()
            if np.any(data < 0):
                print(f"[check] ✗ Layer '{layer}' has negative values — expected non-negative moments.")
                ok = False
            else:
                print(f"[check] ✓ Layer '{layer}': shape {data.shape}, "
                      f"range [{data.min():.2f}, {data.max():.2f}]")

    # Check obs
    if "X_umap" not in adata.obsm:
        print("[check] ⚠ No UMAP found (adata.obsm['X_umap']) — run sc.tl.umap() for visualisation.")

    if ok:
        print("[check] ✓ AnnData is ready for GyreVI.setup_anndata().")

    return ok
