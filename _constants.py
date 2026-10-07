"""Shared constants for GyreVI.

Registry keys and default cell cycle gene lists (human + mouse).

Core CC genes for L_cc are no longer a hand-picked subset — get_cc_genes()
returns the full Seurat G1/S + G2/M union as both cc_genes and core_cc_genes.
The per-gene R_min threshold in the module adapts to each gene's expression
level, so lowly expressed or non-cycling genes are automatically ignored.
"""

# ---------------------------------------------------------------------------
# AnnData registry keys — must match the strings used in setup_anndata()
# and the module's tensor dict access.
# ---------------------------------------------------------------------------
SPLICED_KEY   = "spliced"    # Ms — KNN-smoothed normalised spliced
UNSPLICED_KEY = "unspliced"  # Mu — KNN-smoothed normalised unspliced
BATCH_KEY     = "batch"      # integer batch index


# ---------------------------------------------------------------------------
# Default cell cycle gene lists
#
# Sources:
#   Human — Tirosh et al. 2016 (Seurat cc.genes); Regev lab hallmarks
#   Mouse — converted from human via MGI ortholog mapping
#
# These are used for two purposes:
#   1. Gene gate initialisation (+3.0 bias for these genes)
#   2. L_cc hinge loss (core subset only)
#
# Usage:
#   from _constants import CC_GENES_HUMAN, CORE_CC_GENES_HUMAN
#   from _constants import CC_GENES_MOUSE, CORE_CC_GENES_MOUSE
# ---------------------------------------------------------------------------

# --- Human cycling genes (G1/S + G2/M phases) ---
CC_GENES_HUMAN = [
    # G1/S
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
    # G2/M
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
    "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


# --- Mouse cycling genes (ortholog-mapped from human) ---
CC_GENES_MOUSE = [
    # G1/S
    "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung",
    "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip",
    "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76",
    "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51",
    "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm",
    "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8",
    # G2/M
    "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a",
    "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf",
    "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb",
    "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp",
    "Cdca3", "Hn1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1",
    "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr",
    "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf",
    "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa",
]


def get_cc_genes(organism: str = "human") -> tuple[list[str], list[str]]:
    """Return (cc_genes, core_cc_genes) for the given organism.

    Both returned lists are identical — the full Seurat G1/S + G2/M union.
    Core CC genes are no longer a hand-picked subset: the per-gene R_min
    threshold in GyreVIModule scales each gene's hinge floor by its own
    mean expression, so lowly expressed or non-cycling genes in this dataset
    are automatically ignored (R_min_g ≈ 0 → hinge never activates).

    Parameters
    ----------
    organism : str
        ``'human'`` or ``'mouse'`` (case-insensitive).

    Returns
    -------
    tuple[list[str], list[str]]
        (cc_gene_list, core_cc_gene_list)  — both are the full CC list

    Example
    -------
    >>> cc_genes, core_genes = get_cc_genes("mouse")
    >>> model = GyreVI(adata, cc_gene_names=cc_genes, core_cc_gene_names=core_genes)
    """
    org = organism.lower()
    if org in ("human", "hs", "homo_sapiens"):
        return CC_GENES_HUMAN, CC_GENES_HUMAN
    elif org in ("mouse", "mm", "mus_musculus"):
        return CC_GENES_MOUSE, CC_GENES_MOUSE
    else:
        raise ValueError(
            f"Unknown organism '{organism}'. Choose 'human' or 'mouse'."
        )
