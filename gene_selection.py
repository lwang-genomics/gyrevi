"""Data-driven selection of cycling genes from (s, u) geometry.

Chooses which genes get the constrained (z-bypassing) baseline WITHOUT a curated
cell cycle list, so the list becomes a result to compare against rather than an
input the model depends on.

WHAT WAS EXPECTED, AND WHY IT WAS WRONG
---------------------------------------
The design notes argue cycling genes trace closed *loops* in (s, u) while
differentiation genes trace open arcs, which predicts cycling genes should show a
large u/s phase lag (low corr(s,u) = cos(delta)) and two u-branches at a given s.

Measured on the fibroblast benchmark, both statistics came out ANTI-predictive
(AUROC 0.175 and 0.216 against the Tirosh list). The premise is inverted: since
delta = arctan(omega / gamma_g) and degradation is fast relative to the cycle for
most genes, delta is near 0, so the loop is extremely thin — cycling genes look
like tight lines, not fat loops. Loop AREA scales with sin(delta) and is simply
not a usable signal here.

WHAT ACTUALLY DISCRIMINATES
---------------------------
The inverse: **tight s-u coupling**. Cycling genes are coherently driven, so their
spliced and unspliced levels track each other cleanly (|corr| 0.693 vs 0.382 for
the rest, AUROC 0.825). This is not an expression artifact — mean expression alone
gives AUROC 0.517. Ring-hollowness adds a little on its own (0.645).

So the selector measures COHERENT KINETICS, not loop closure. Named accordingly.
"""

from __future__ import annotations

import numpy as np

EPS = 1e-8


def _as_dense(X):
    return np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float64)


def _zscore(X):
    return (X - X.mean(0)) / (X.std(0) + EPS)


def coupling_score(Ms, Mu):
    """|corr(s, u)| per gene. High = coherent kinetics = cycling-like.

    AUROC 0.825 against the Tirosh list on the fibroblast benchmark. Note this is
    the INVERSE of what the loop-closure argument predicts; see the module
    docstring for why that prediction fails.
    """
    s, u = _zscore(_as_dense(Ms)), _zscore(_as_dense(Mu))
    return np.abs((s * u).mean(0))


def hollowness_score(Ms, Mu):
    """Ring-likeness of the (s, u) cloud: low central density relative to the rim.

    Weakly informative on its own (AUROC 0.645) — kept because it is the one
    genuinely loop-shaped signal that survived, and it is not redundant with
    coupling.
    """
    s, u = _zscore(_as_dense(Ms)), _zscore(_as_dense(Mu))
    r = np.sqrt(s ** 2 + u ** 2)
    return np.quantile(r, 0.1, axis=0) / (np.quantile(r, 0.9, axis=0) + EPS)


def branch_score(Ms, Mu, n_bins: int = 20):
    """Within-s-bin variance of u, as a fraction of u's total variance.

    Cells are binned by their rank in s, so this asks: at a given spliced level,
    how much does unspliced still vary? A closed loop has two limbs and so a
    large within-bin spread; a single-valued arc has almost none.
    """
    s, u = _zscore(_as_dense(Ms)), _zscore(_as_dense(Mu))
    n, g = s.shape
    order = np.argsort(s, axis=0)
    u_sorted = np.take_along_axis(u, order, axis=0)
    edges = np.linspace(0, n, n_bins + 1).astype(int)
    within = np.zeros(g)
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < 3:
            continue
        within += u_sorted[a:b].var(0) * (b - a)
    return within / n            # u is z-scored, so total variance is 1


def expression_filter(Ms, min_mean: float = 0.1, min_cv: float = 0.1):
    """Genes with enough signal for the geometry to mean anything."""
    X = _as_dense(Ms)
    mean = X.mean(0)
    cv = X.std(0) / (mean + EPS)
    return (mean > min_mean) & (cv > min_cv)


def cycling_scores(Ms, Mu, n_bins: int = 20):
    """All statistics plus a combination, with the expression filter applied."""
    keep = expression_filter(Ms)
    coupling = coupling_score(Ms, Mu)
    hollow = hollowness_score(Ms, Mu)
    branch = branch_score(Ms, Mu, n_bins=n_bins)

    def rank01(v):
        r = np.empty_like(v)
        r[np.argsort(v)] = np.arange(v.size)
        return r / max(v.size - 1, 1)

    # 0.9/0.1, not 0.75/0.25 — hollowness is PHASE-BIASED and at 25% it silently
    # excluded an entire cell cycle phase.
    #
    # hollowness was intended as a ring-shape statistic, but it correlates +0.304 with the
    # fraction of cells above a gene's own mean: it rewards genes that are ON for a large
    # fraction of the cycle. Mitosis is brief, so G2/M genes sit near the origin in most
    # cells, look "filled in", and are demoted. Measured on Battich, top-100 by each score:
    #
    #   statistic                     AUROC   Tirosh   G2/M markers   S markers
    #   coupling only                 0.848   21/100        5/9           1/9
    #   0.75 coupling / 0.25 hollow   0.852   23/100        0/9           3/9   <- was default
    #   0.90 coupling / 0.10 hollow   0.854   28/100        3/9           2/9   <- now default
    #
    # The old weighting bought +0.004 AUROC and cost every canonical G2/M marker. TOP2A has
    # the highest coupling of any gene tested (0.897) and still fell to rank 318. Since the
    # discovery claim rests on this selector, a score that cannot see M phase is a real
    # defect, not a tuning preference.
    combined = 0.9 * rank01(coupling) + 0.1 * rank01(hollow)
    out = {"coupling": coupling, "hollowness": hollow, "branch": branch,
           "combined": combined, "expressed": keep}
    for k in ("coupling", "hollowness", "branch", "combined"):
        out[k] = np.where(keep, out[k], 0.0)
    return out


def select_cycling_genes(Ms, Mu, n_top: int = 100, statistic: str = "coupling",
                         n_bins: int = 20):
    """Indices of the top-n genes by the chosen statistic."""
    sc = cycling_scores(Ms, Mu, n_bins=n_bins)[statistic]
    return np.argsort(sc)[::-1][:n_top]
