# GyreVI

**Unsupervised recovery of a latent circular phase from sparse, high-dimensional count data.**

A variational autoencoder that assigns every observation an angle on a circle — a continuous latent
phase — from a ~2,000-dimensional sparse count vector, with no labels, and separates that phase from
a residual state vector. The phase is validated against an independent physical measurement of the
same quantity.

The application is cell-cycle inference in single-cell RNA sequencing. The statistical problem is
general: **recover a periodic latent state from noisy high-dimensional observations, when the only
ground truth is itself a noisy proxy and every headline metric has a trivial estimator that has to be
ruled out first.**

Work done at Institut Curie. Manuscript in preparation — see
[Release status](#release-status) for what is and is not in this repository.

---

## The problem

| | |
|---|---|
| Observations | 5,422 cells × ~2,000 genes, sparse integer counts |
| Latent target | a phase θ ∈ [0, 2π) per observation, plus a residual state vector **z** |
| Supervision | none — the phase is recovered unsupervised |
| Ground truth | an independent two-channel fluorescent protein readout, itself noisy |

Three things make this harder than the dimensionality suggests.

**The target is circular.** Means, variances, correlation and regression are all undefined or
misleading on angles. 0.01 and 6.27 radians are neighbours. Every comparison has to go through
circular statistics, and the estimator has to be invariant to the arbitrary rotation and handedness
of the inferred frame — which means the alignment step itself can manufacture a result if it is not
constrained.

**Ground truth is a proxy, not a label.** The fluorescent readout measures two proteins whose
abundance tracks the phase; it does not measure the phase. It has a window where both channels are
dark and the angle is genuinely undefined, and its own circular spread (314°) is a ceiling on what
agreement can mean.

**Every obvious metric has a trivial baseline.** Two of the five parameters this model reports can be
approximated to within a few percent by an estimator with no model in it at all. Those two are
reported *net* of that baseline, and the two with no trivial baseline are the ones the method claim
actually rests on. See [Validation](#validation-and-negative-results).

---

## Results

Six seeds, paired, measured against the protein ground truth on the same 5,422 cells.

| quantity | measured | trivial baseline | verdict |
|---|---|---|---|
| **δ** per-gene lag between the two count channels | **+0.598 ± 0.029** | *none exists* | **the load-bearing claim** |
| **φ** per-gene peak phase | **+0.706 ± 0.044** | *none exists* | holds |
| **λ** per-gene channel ratio | +0.980 ± 0.000 | **+0.965** (ratio of means) | holds, but nearly free |
| **R** per-gene amplitude | +0.912 | **+0.828** (mean expression) | +0.749 with expression partialled out |
| **θ** per-cell phase | +0.473 ± 0.009 | **+0.463** (two gene-set means) | **treat as a tie** |

### Phase recovery against published methods

| | method | Spearman vs ground truth | median angular error |
|---|---|---|---|
| *oracle* | scVI latent + ridge probe **fitted to the labels** | *+0.669* | *15°* |
| label-free | two gene-set means + `atan2` (Seurat) | +0.463 | 31° |
| label-free | **GyreVI θ** | **+0.473 ± 0.009** | 31° |
| label-free | CycleVI | +0.438 / +0.573 | 28° |
| label-free | DeepCycle | +0.275 | 49° |

Three honest qualifications, because the table does not speak for itself:

1. **The top row is a ceiling, not a competitor.** It is fitted to the ground-truth labels. Nothing
   label-free should be expected to reach it, and its presence is the point: it bounds how much phase
   information the data contains at all.
2. **θ does not beat the simple baseline.** +0.473 ± 0.009 against +0.463 is inside run-to-run drift.
   It is reported here as a tie, and the method claim rests on δ and φ instead — quantities the
   baseline cannot produce at all.
3. **CycleVI has two numbers because the metric convention decides the ranking.** +0.438 under
   circular alignment, +0.573 under a best-shift search; the authors publish 0.513. These are the
   same predictions scored three ways. A single number here would be a choice presented as a result,
   so both are given.

---

## Validation and negative results

This is the part of the project worth reading. Every item below silently produced a wrong number
before it was caught, and each is recorded so it is not retried. Most are not specific to biology.

**Selection and filtering**

- **A confidence filter needs a matched baseline.** Restricting to the model's most confident half
  lifts θ from +0.557 to +0.746 — but lifts the *simple baseline* from +0.463 to +0.748 on the same
  cells. The filter selects easy observations; it does not demonstrate a better estimator.
- **A subset correlation is only meaningful if the subset still spans the circle.** Filtering by
  signal brightness looks like it improves agreement, but collapses the retained arc from ~314° to
  68°. The correlation rises because the range shrank. Filtering by model confidence keeps ~300° and
  is safe. *Restriction of range was the single most repeated trap in this project.*
- **Never subset by a fitted quantity to rescue a claim.** Selecting on the model's own output
  inflated one parameter from +0.278 to +0.476. Model-independent selection is not circular;
  model-dependent selection is.
- **Held-out AUROC is not comparable across selectors.** A better selector pulls easy positives into
  its own selected set and is then scored on a harder remainder. Score both on the identical mask —
  and on the identical matrix: the same predictions score 0.933 on raw counts and 0.908 on smoothed
  counts, because smoothing lifts every fit.

**Controls and baselines**

- **A control that correlates with the thing it controls for is not a control.** One claim was
  retracted on a table whose two columns turned out to be the same finding printed twice. The
  agreement read as corroboration and was arithmetic.
- **Check the units of every baseline.** One published comparator reports phase on [0, 1] and another
  in radians. Passing the first straight into a radian comparison scored it +0.100 over a 49° arc
  instead of +0.275 over 331°. The absurd arc was the tell, and it sat printed in a results table for
  a full session before anyone read it.
- **Never tune on a validation set you have not checked the sign of.** A second dataset used for
  tuning turned out to be *anti*-correlated with truth: the configuration that scored best on the
  validated data was the worst on it. Roughly 14 runs had been tuned against it.

**Inference and comparison**

- **Only paired comparisons support inference here.** The hardware backend is non-deterministic, so
  even fixed seeds vary. One conclusion in this project was drawn from three separately-run
  benchmarks; a paired re-test showed both of its claims were false. A single-seed +0.540 became
  +0.453 over five seeds.
- **Averaging only cancels independent errors.** Ensembling six seeds made the phase *worse*
  (+0.457 against the mean single seed's +0.466), because the error is systematic — a misaligned
  objective — not stochastic. A control confirmed the alignment step was not at fault.
- **Prefer an estimator that cannot be wrong about its own convergence.** Three wrong numbers came
  from reading a model's own fitted output as if it were a measurement. A closed form disagreed with
  all three.
- **Verify a new code path is reached before trusting a null result.** A failed patch left two arms
  running identical code; they came out near-identical and that was reported as a refutation.
  Near-identical arms are a symptom, not a result.
- **A preprocessing step can degrade the thing it normalises.** A rank-uniformisation step *repairs*
  one badly-distributed input (+0.426 → +0.463) and *destroys* an already-even one (+0.579 → +0.485),
  because the true density is genuinely non-uniform. Check a preprocessing step against the input it
  will actually receive, not the one it was written for.

**A bug that passed every health check**

The reference used to orient the circular frame was reversed in handedness. It survived the entire
test suite, because every check either sign-corrects internally or compares two quantities that flip
together. It was caught only by comparing against an external reference that could not flip. Fixing
it forced the retraction of a reported gain *and* a reported cost — two artifacts of the same bug,
pointing in opposite directions.

---

## What is in this repository

```
circular.py                  circular statistics and the phase-comparison metrics
                             (wrap, circular mean, resultant length, rotation-invariant
                             alignment, circular correlation) — standalone, no dependencies
                             on the model
preprocess.py                count preprocessing
prepare_battich.py           dataset preparation and the four measurement layers
gene_selection.py            model-independent gene selection
baselines/scvi_baseline.py   the supervised-probe ceiling
baselines/run_cyclevi.py     published comparator, re-run locally under matched conditions
figures/                     benchmark results
```

`circular.py` is the piece most likely to be useful on its own: the alignment and scoring routines
are what make the comparisons in this repository rotation- and handedness-invariant, and they are
where several of the traps above were eventually caught.

## Release status

The model implementation (`gyrevi_module.py`, `gyrevi_model.py`) and the architecture notes are
**held back pending the manuscript**, along with the development history. This repository contains
the evaluation harness, the baselines and the results — deliberately, since that is the part that is
reusable and the part the claims rest on.

An earlier version of this work was independently published by another group during development, so
the method is being released with the preprint rather than before it. The full implementation and
commit history are available on request.

## Data

The dataset is Battich et al. RPE1-FUCCI (published, publicly available). Raw data is not vendored
here; `prepare_battich.py` documents the retrieval and preparation. Comparator outputs under
`data/` are reproductions of published methods, with their provenance recorded alongside.

## License

MIT — see [LICENSE](LICENSE).
