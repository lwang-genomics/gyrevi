# DeepCycle baseline on Battich RPE1-FUCCI

`deepcycle_theta.csv` — DeepCycle's `cell_cycle_theta` on [0,1] for all 5,422 cells, so the
head-to-head does not need to be re-run.

## Reproducing

DeepCycle pins Python 3.7 / TensorFlow 2.2.0-gpu / CUDA 10.1 (`DeepCycle_env.yml`), none of which
has an arm64 macOS build. It runs unmodified on **TF 2.15** — the last release defaulting to Keras 2,
which is what this code targets — with one edit (their source is not redistributed here — apply it to a fresh clone):

* `import tensorflow_datasets as tfds` is commented out. It is imported and **never used**, and it
  drags in `tensorflow_metadata`, which needs a protobuf newer than TF 2.15 accepts.

```
uv venv --python 3.11 dcenv
VIRTUAL_ENV=dcenv uv pip install "tensorflow==2.15.1" "tensorflow-probability==0.23.0" \
    "numpy<2" anndata scanpy seaborn matplotlib scipy scikit-learn
dcenv/bin/python DeepCycle.py --input_adata   # after commenting out the tfds import rpe1_gyrevi.h5ad \
    --gene_list gene_list_present.txt --base_gene TOP2A \
    --expression_threshold 0.5 --output_adata out_deepcycle.h5ad
```

`gene_list_present.txt` is the 394 genes of DeepCycle's human GO cell_cycle annotation that survive
our HVG subset. `TOP2A` was chosen as `--base_gene` for having the highest CV among well-expressed
canonical cycling genes (0.73 at mean Ms 2.28), which is what its bimodality requirement asks for.

**Expect a traceback at the very end.** DeepCycle's optional "PRELIMINARY ANALYSIS" plot reads
`adata.obsm['X_umap']`, which our AnnData does not carry. It raises **after** the output h5ad is
written, so the result is complete — exit code is 0 and `cell_cycle_theta` is present.

Training took ~16 min on CPU (best epoch 152, MSE 0.669).
