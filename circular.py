"""Circular statistics for comparing cell cycle phase estimates.

Phase estimates from different methods are only defined up to two nuisance
parameters:

1. **Rotation** — there is no canonical zero point on the cycle, so every
   method picks an arbitrary origin. (Mozdzanowski et al. handle this by
   circularly shifting DeepCycle by 2.9 rad before plotting.)
2. **Direction** — nothing in an unsupervised model fixes which way round the
   circle time runs, so one method's theta may be the other's -theta.

Any comparison must therefore align on (direction, rotation) *first* and only
then measure agreement. Correlating ``cos(theta)`` directly, a common shortcut,
is wrong on both counts: ``cos`` is an even function so it cannot tell theta
from -theta, and it silently assumes both methods share an origin.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr

TWO_PI = 2.0 * np.pi


# ── basic helpers ────────────────────────────────────────────────────────────

def wrap_pi(theta):
    """Wrap angles to (-pi, pi]."""
    return (np.asarray(theta, dtype=float) + np.pi) % TWO_PI - np.pi


def wrap_2pi(theta):
    """Wrap angles to [0, 2pi)."""
    return np.asarray(theta, dtype=float) % TWO_PI


def to_radians(x, period=1.0):
    """Rescale a phase expressed on ``[0, period)`` to radians on ``[0, 2pi)``.

    DeepCycle reports its phase on [0, 1], not in radians — comparing it
    against a radian-valued theta without this conversion silently compares
    two different units.
    """
    return wrap_2pi(np.asarray(x, dtype=float) * (TWO_PI / period))


def circ_mean(theta):
    """Circular mean of angles (radians)."""
    theta = np.asarray(theta, dtype=float)
    return np.angle(np.mean(np.exp(1j * theta)))


def resultant_length(theta):
    """Mean resultant length R in [0, 1] — 1 = perfectly concentrated."""
    theta = np.asarray(theta, dtype=float)
    return float(np.abs(np.mean(np.exp(1j * theta))))


def quantile_uniform_angle(values):
    """Rank-transform values to angles uniform on [0, 2pi).

    This is the transform Mozdzanowski et al. apply to the Seurat-derived
    angle so it covers the circle evenly before it is used as a prior. Ties
    are broken by average rank.
    """
    values = np.asarray(values, dtype=float)
    n = values.size
    order = values.argsort(kind="mergesort")
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(n, dtype=float)
    return wrap_2pi(ranks / n * TWO_PI)


def seurat_angle(s_score, g2m_score, quantile=True):
    """Seurat reference angle: ``atan2(G2M, S)``, optionally rank-uniformised.

    Mirrors the reference-angle construction in Mozdzanowski et al. §2.1.3.
    """
    raw = np.arctan2(np.asarray(g2m_score, dtype=float),
                     np.asarray(s_score, dtype=float))
    return quantile_uniform_angle(raw) if quantile else wrap_2pi(raw)


# ── alignment ────────────────────────────────────────────────────────────────

def align_circular(reference, query):
    """Align ``query`` onto ``reference`` over direction and rotation.

    For each direction ``s`` in {+1, -1} the rotation that best matches the
    reference has a closed form: the circular mean of ``reference - s*query``.
    We take whichever direction leaves the tighter residual.

    Returns
    -------
    aligned : ndarray
        ``wrap_2pi(s * query + shift)``.
    info : dict
        ``sign``, ``shift`` (radians), and ``R`` — the mean resultant length of
        the residual, in [0, 1], which is the agreement after alignment.
    """
    reference = np.asarray(reference, dtype=float)
    query = np.asarray(query, dtype=float)

    best = None
    for sign in (1.0, -1.0):
        resid = reference - sign * query
        shift = circ_mean(resid)
        R = resultant_length(resid - shift)
        if best is None or R > best["R"]:
            best = {"sign": sign, "shift": float(shift), "R": float(R)}

    aligned = wrap_2pi(best["sign"] * query + best["shift"])
    return aligned, best


# ── association measures ─────────────────────────────────────────────────────

def circ_corr(a, b):
    """Jammalamadaka-Sarma circular correlation coefficient, in [-1, 1].

    Rotation-invariant by construction (it centres both inputs on their
    circular means) and O(n), unlike the O(n^2) Fisher-Lee coefficient.
    Sign-sensitive: a reversed cycle gives a negative value.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    sa = np.sin(a - circ_mean(a))
    sb = np.sin(b - circ_mean(b))
    denom = np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2))
    return float(np.sum(sa * sb) / denom) if denom > 0 else np.nan


def compare_phases(reference, query, ref_name="reference", query_name="query"):
    """Full comparison of two phase estimates.

    Aligns ``query`` onto ``reference`` first, then reports agreement. The
    ``spearman`` field is computed on the aligned angles and so is the value
    directly comparable to the Spearman correlations quoted in the CycleVI
    preprint.
    """
    reference = wrap_2pi(np.asarray(reference, dtype=float))
    query = wrap_2pi(np.asarray(query, dtype=float))

    ok = np.isfinite(reference) & np.isfinite(query)
    reference, query = reference[ok], query[ok]

    aligned, info = align_circular(reference, query)
    rho, pval = spearmanr(reference, aligned)

    # Median absolute angular error, in degrees — an interpretable companion
    # to R that says how far off a typical cell is.
    err = np.abs(wrap_pi(reference - aligned))

    # Degeneracy guard.
    #
    # R alone cannot tell a real match from a collapsed prediction. If the query
    # is (near-)constant, the residual is just `reference - const`, so R comes
    # out equal to the REFERENCE's own concentration no matter what the model
    # learned. On the fibroblast benchmark a literal constant scored R=0.510
    # against DeepCycle, versus 0.514 for a genuinely trained-but-collapsed θ —
    # indistinguishable. `R_null` is that constant-prediction score, and
    # `R_excess = R - R_null` is the part actually attributable to the model.
    ref_conc = resultant_length(reference)
    qry_conc = resultant_length(query)
    R_null = ref_conc
    degenerate = bool(qry_conc > 0.95)

    return {
        "reference": ref_name,
        "query": query_name,
        "n": int(ok.sum()),
        "circ_corr": circ_corr(reference, aligned),
        "spearman": float(rho),
        "spearman_p": float(pval),
        "R": info["R"],
        "R_null": float(R_null),
        "R_excess": float(info["R"] - R_null),
        "ref_concentration": float(ref_conc),
        "query_concentration": float(qry_conc),
        "degenerate": degenerate,
        "median_abs_err_deg": float(np.degrees(np.median(err))),
        "direction": "same" if info["sign"] > 0 else "reversed",
        "shift_rad": info["shift"],
        "aligned": aligned,
        "mask": ok,
    }
