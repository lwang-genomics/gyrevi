# GyreVI

**A variational autoencoder that jointly infers cyclical dynamics and cellular state from single-cell RNA-seq data.**

From ~2,000 sparse gene counts per cell, GyreVI places each of 5,422 cells on a circle, its position in the cell cycle. For every gene, it estimates when the gene peaks and how far
its unspliced transcripts run ahead of the spliced ones. Every output is scored against an independent
protein measurement of each cell's position.

![The data: protein ground truth, and a cycling gene's spliced and unspliced RNA](figures/data_intro.png)

**a** Two fluorescent reporters give each cell's true position in the cycle; they are never a model input, only the ground truth for scoring.
**b–c** In NUF2, one of the clearest cycling genes, new (unspliced) RNA rises before mature (spliced) RNA, so the two trace a loop.
The peak φ and the lead δ (25° here, about 10° for a typical cycling gene) are what GyreVI estimates for every gene.

![The model at a glance](figures/model_schematic.png)

The encoder maps each cell's counts to θ, its position on the circle, and z, the rest of its state. The
decoder turns them back into expected counts, with a peak φ, a lead δ and an amplitude R for each gene.

**Result: a simple linear baseline sets the bar, and shows how much signal is left to capture.** A phase
from PCA on 90 known cell-cycle genes, then a least-squares cosine fit per gene, leads on every quantity.

![GyreVI's learned estimates against simple estimators](figures/simple_vs_model.png)

| correlation with ground truth | GyreVI, learned | least squares at GyreVI's phase | CycleVI phase + least squares | Seurat phase + least squares | PCA phase + least squares |
|---|---|---|---|---|---|
| cell phase θ | +0.47 ± 0.01 | — | +0.44 | +0.46 | **+0.58** |
| gene peak phase φ | +0.71 ± 0.04 | +0.77 ± 0.01 | +0.77 | +0.77 | **+0.82** |
| unspliced lead δ | +0.60 ± 0.03 | +0.82 ± 0.02 | +0.93 | +0.93 | **+0.95** |
| gene amplitude R | +0.91 ± 0.01 | +0.97 ± <0.01 | +0.97 | +0.96 | **+0.98** |

Spearman correlation (θ after circular alignment), except φ (circular correlation). GyreVI: mean ± sd over 6
paired seeds. The last three columns reproduce with `python baselines/two_step.py`; the CycleVI column
needs CycleVI's output from `baselines/run_cyclevi.py` first.

The benchmark pinpoints two gaps:

- **Phase.** GyreVI is trained to follow the PCA phase as a reference, and ends up further from the truth
  than that reference (+0.47 vs +0.58).
- **Gene parameters.** Even at the model's own phase, least squares beats its learned parameters
  (δ: +0.82 vs +0.60).

**Why this is useful.** The PCA result shows the data hold more phase information than the model extracts
(+0.58; a supervised probe reaches +0.68), so the gap is in the model, not the data. The two gaps above
say where to look. A natural next step is to anchor the model at the PCA solution and learn only
corrections.

*Caveat:* the per-gene targets are themselves cosine fits at the protein phase, on the same counts, which
favours estimators of that form. Kinetics measured independently by metabolic labelling, which this
dataset supports, would be the cleaner test.

## Phase against published methods

![Cell-phase accuracy of each method](figures/phase_methods.png)

Circular alignment for every method. The oracle is fitted to the labels: a ceiling, not a competitor.
GyreVI and the oracle: mean of 6 seeds. With the best of 144 rotations instead (chosen using the labels),
every method scores higher and the middle three change order; the PCA leads under both.

![Inferred against protein phase, per cell](figures/phase_scatter.png)

Top row: deep models; bottom: simple baselines. GyreVI: the seed closest to the 6-seed mean. The corners
are the same point on the circle (0 = 2π).

## Confidence

![Calibration of each confidence score, and how much of the cycle filtering keeps](figures/confidence.png)

**a** GyreVI's confidence ρ is calibrated: agreement rises steadily with it. The PCA's radius behaves the
same way, and sets the bar here too. **b** Filtering on either keeps most of the cycle. Filtering on signal
brightness keeps a 68° sliver, and agreement falls to zero, because brightness tracks the phase itself.

## Checks behind every number

- A simple estimator for every learned quantity, to set the target to beat.
- Paired seeds only: the GPU backend is non-deterministic.
- Range before correlation: a filter correlated with the target truncates it (panel b).
- An external reference for every internal check: a reversed phase frame once passed all of them.

## Repository

| | |
|---|---|
| `baselines/two_step.py` | the simple estimators above, from this repository alone |
| `circular.py` | circular statistics: alignment, circular correlation, resultant length |
| `prepare_battich.py` · `preprocess.py` · `gene_selection.py` · `_constants.py` | data preparation |
| `baselines/scvi_baseline.py` · `baselines/run_cyclevi.py` | the oracle, and CycleVI under matched conditions |

GyreVI's own code is not yet public; its numbers come from a 6-seed benchmark. Work done at Institut Curie.
Data: Battich et al., RPE1-FUCCI (public). License: MIT.
