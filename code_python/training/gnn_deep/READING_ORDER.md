# polykernels: models and pipeline

This is a snapshot (26 October 2026) of the code behind the regression results we exported.
It covers four things:
- how molecules and tabular inputs become model inputs;
- the six models;
- how they are trained and evaluated;
- the three attribution methods.

The files are listed below in a suggested reading order. The other files in the folder are
support code that these import.

## Overview

Each row is a molecule (SMILES) plus tabular inputs, the dataset's `x_cols`. These are
process conditions and material or solvent descriptors.
- **Graph models** see the molecule as a graph, with the `x_cols` attached as graph-level
  globals.
- **Tabular models** see only the `x_cols`.

The results come from two kinds of run:
- **IID:** 3 seeds × 5-fold cross-validation.
- **OOD:** leave-one-cluster-out splits.

The attributions explain each IID fold's model on that fold's own held-out rows.

## The six models

| Code name | Export name | What it is |
|---|---|---|
| `gnn_mlp` | GNN-MLP | GNN encoder and MLP head; point prediction |
| `gnn_gp` | GNN-GP | The same GNN, with a sparse variational GP head per target on its embedding |
| `gnn_gpsn` | GNN-GP-SN | GNN-GP with a spectral-norm-bounded, distance-preserving encoder (SNGP / DUE style) |
| `mlp_gp` | MLP-GP | MLP over the `x_cols` only, with a GP head |
| `mlp_gpsn` | MLP-GP-SN | Spectral-norm residual MLP, with a GP head |
| `gnn_mlpens` | GNN-MLP-Ens | Deep ensemble of 5 GNN-MLPs, each predicting a mean and a variance |

## Reading order

### 1. Data → model inputs
1. [polykernels/regression_data.py](polykernels/regression_data.py): the seven
   regression problems, with their targets, `x_cols` and cluster columns.
2. [polykernels/molgraph/templates.py](polykernels/molgraph/templates.py), then
   [mol_to_graph.py](polykernels/molgraph/mol_to_graph.py): how SMILES become graphs.
   - The `CUSTOM` template is the one used.
   - [features.py](polykernels/molgraph/features.py) defines the atom and bond features.
     It is long, so skim it.
3. [polykernels/datasets.py](polykernels/datasets.py): the graph dataset.
   `overwrite_dataset_globals` writes the standardized `x_cols` into each graph's globals
   `u`.

### 2. Models
4. [models/mlp.py](polykernels/models/mlp.py): the MLPs, used inside the GNN and as the
   tabular feature extractor.
5. [models/gnn.py](polykernels/models/gnn.py): the GNN, in three stages:
   - an encoder;
   - weight-tied graph-network layers: an edge update, GATv2 node attention, and a pooled
     global update. The globals condition every layer through AdaLN.
   - a readout.
6. [models/spectral.py](polykernels/models/spectral.py): the spectral-norm bound used by
   the `*_gpsn` models.
7. [models/gp.py](polykernels/models/gp.py): the GP heads, how their inducing points are
   seeded, and ELBO training.
8. [models/ensemble.py](polykernels/models/ensemble.py): GNN-MLP-Ens.
9. [models/losses.py](polykernels/models/losses.py) and
   [models/training.py](polykernels/models/training.py): the losses and training loops.

### 3. Training and evaluation
10. [scripts/train_gnn_regression.py](scripts/train_gnn_regression.py): the entry point.
    - `build_model` builds each model type from its hyperparameters.
    - `warmup_lazy_modules` does the GP-specific setup before training.
    - `train_one_fold` and `evaluate_fold` train and score one fold.
    - `make_cv_jobs` lays out the IID and OOD jobs.
11. [polykernels/splits.py](polykernels/splits.py): k-fold, leave-one-cluster-out and
    size-matched IID splits.
12. [polykernels/uq_metrics.py](polykernels/uq_metrics.py): RMSE, MAE, R², NLL, and the
    calibration metrics.

### 4. Attribution
13. [scripts/analyze_shap.py](scripts/analyze_shap.py): SHAP over the `x_cols`, for any
    model. For graph models, each row's molecule is held fixed.
14. [polykernels/gpshap.py](polykernels/gpshap.py), then
    [scripts/analyze_gpshap.py](scripts/analyze_gpshap.py): GP-SHAP, which gives
    stochastic Shapley values with credible intervals, for the GP models.
15. [polykernels/attribution.py](polykernels/attribution.py), then
    [scripts/analyze_attributions.py](scripts/analyze_attributions.py): Integrated
    Gradients for the GNNs, over atoms and node features and over the `x_cols` (globals).
16. [scripts/attribution_consistency.py](scripts/attribution_consistency.py) and
    [scripts/attribution_agreement.py](scripts/attribution_agreement.py): Kendall's W
    of the importance rankings across the 15 CV fold models (5 folds × 3 seeds), and
    between methods. Per-row and per-atom W compare the 3 seeds instead, since each row is
    held out by one fold model per seed.

## Support files
These are imported by the files above but aren't needed to follow the models:
- `polykernels/__init__.py`
- `utils.py`, which includes `set_seed`: seeding plus deterministic GPU kernels
- `paths_and_constants.py`, which sets the file locations
- `vis.py`, `metrics.py`, `chem.py`, `io.py`, `stratification.py` and `types.py`
- `molgraph/__init__.py` and `molgraph/graphwave.py`

## Running it
None of this is needed for reading. The data, the tuned hyperparameters and the cluster job
scripts are not included.

- **Environment:** create it from `environment.yml`, then run `pip install -e .` in this
  folder. The results used RDKit 2026.03.4; other RDKit versions can change some bond
  features, and so the graph inputs.
- **Data:** the dataset CSVs go in `data/processed/`, under the file names in
  `polykernels/paths_and_constants.py`.
- **Hyperparameters:** the results used hyperparameters tuned for each problem and model.
  Without `--hparams_path`, the script uses the defaults in `get_default_hparams`.
- **A quick test run:**
  ```
  python scripts/train_gnn_regression.py --problem dataset5_rg --model_type gnn_gp \
      --split_types '[iid]' --seeds '[6]' --n_folds 2 --results_dir results/test
  ```
- **OOD splits:** `--ood_clusters structure_cluster` (or `continuous_cluster`) needs those
  cluster labels imported next to the data.
- **Attribution:** the attribution scripts explain a saved model. They expect
  `model.safetensors`, `hparams.yaml` and `standardization.npz` in
  `<results_dir>/<problem>_<model_type>/`. For the reported results, those were the
  fold models saved by `train_gnn_regression.py --save_models true`.
