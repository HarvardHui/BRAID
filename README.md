# BRAID and BEGONE

**BRAID** (Batch-effect Reconciliation and Anomaly-aware Imputation for Data) is an evaluation framework that measures how batch effects change the performance of missing value imputation methods. It compares each method under a batch-free condition (batch-ve) and a batch-affected condition (batch+ve) with identical missing value positions, scoring both against the same ground truth.

**BEGONE** (Batch Effect Gradient Observation on Null Estimation) extends BRAID by injecting simulated batch effects over a grid of additive and multiplicative severities. It shows how each method degrades as batch effects grow stronger.

This repository contains the R code for both frameworks and the scripts that reproduce every figure and table in the manuscript:

> Hui HWH, Goh WWB. BRAID: Batch-effect Reconciliation and Anomaly-aware Imputation for Data. *Bioinformatics Advances*.

## Repository layout

```
BRAID/
├── R/BRAID_functions.R          function library, sourced by every script
├── main_figures.R               every main-manuscript figure, in one script
├── run_BRAID.R                  one BRAID setting on its own (optional)
├── run_BEGONE.R                 one BEGONE sweep on its own (optional)
├── supplementary/
│   ├── S01_batch_effect_strength.R   Figures S1–S2
│   ├── S02_feature_subsampling.R     Figure S3
│   ├── S03_mice_rf_limitations.R     Figures S4–S6
│   ├── S04_braid.R                   Figures S7–S20
│   └── S05_begone.R                  Figures S21–S28, Tables S1–S2
├── data/                        input data (see "Input data")
└── results/                     created on first run
```

All scripts find the repository root themselves, so they can be run from the repository folder, from any other folder with `Rscript /path/to/BRAID/run_BRAID.R`, or from RStudio after opening `BRAID.Rproj`. No paths need to be edited.

## Requirements

R (≥ 4.3) and the following packages.

| Source | Packages |
|---|---|
| Bioconductor | `sva`, `limma`, `impute`, `pcaMethods`, `SummarizedExperiment` |
| CRAN | `imputeLCMD`, `mice`, `missForest`, `ggplot2`, `ggsignif`, `car`, `sandwich` |
| Other | `NMFBatch` (Anwar et al., 2026), from https://codeberg.org/AliYoussef/NMFBatch |

```r
install.packages(c("BiocManager", "mice", "missForest", "ggplot2", "ggsignif", "car", "sandwich"))
BiocManager::install(c("sva", "limma", "impute", "pcaMethods", "SummarizedExperiment"))
install.packages("imputeLCMD")
remotes::install_git("https://codeberg.org/AliYoussef/NMFBatch")   # needs the remotes package
```

NMFBatch is optional. If it is not installed, it is skipped with a note and all other methods still run. It is required to reproduce the manuscript results.

The manuscript used `sva` 3.50.0 and `limma` 3.66.0.

## Input data

Place the two datasets in `data/` with exactly these folder and file names:

```
data/
├── Van Puyvelde dataset/
│   ├── HYE5600735_LFQ_FragPipe_design.tsv
│   ├── HYE5600735_LFQ_FragPipe_pro_intensity.tsv
│   ├── HYE6600735_LFQ_FragPipe_design.tsv
│   ├── HYE6600735_LFQ_FragPipe_pro_intensity.tsv
│   ├── HYEqe735_LFQ_FragPipe_design.tsv
│   ├── HYEqe735_LFQ_FragPipe_pro_intensity.tsv
│   ├── HYEtims735_LFQ_FragPipe_design.tsv
│   └── HYEtims735_LFQ_FragPipe_pro_intensity.tsv
└── Haslett & Pescatori datasets/
    ├── DMD-HaslettData.csv
    └── DMD-PescatoriData.csv
```

To keep the data elsewhere, set `BRAID_DATA` to the folder that contains the two dataset folders, e.g. `Sys.setenv(BRAID_DATA = "/path/to/data")` or `export BRAID_DATA=/path/to/data`.

**Van Puyvelde (VP).** FragPipe LFQ output for the DDA runs of a human/yeast/*E. coli* hybrid proteome benchmark, acquired on four instruments (SCIEX TripleTOF 5600, SCIEX TripleTOF 6600+, Orbitrap QE HF-X, Bruker timsTOF Pro). ProteomeXchange identifier PXD028735. Each instrument is one batch.
- `*_design.tsv`: tab-separated, one row per sample, with a `condition` column (A or B). Rows must be in the same order as the sample columns of the matching intensity file.
- `*_pro_intensity.tsv`: tab-separated. Column 1 is `Protein`; column 2 is not used; the remaining columns are sample intensities on the raw (non-log) scale. Protein names must contain `HUMAN`, `YEAST` or `ECOLI`, which define the differential-expression truth (yeast and *E. coli* proteins differ between conditions, human proteins do not).
- Only proteins observed with a non-zero intensity in every run are kept, then log2-transformed.

**Haslett & Pescatori (HP).** Two Duchenne muscular dystrophy microarray studies (Haslett et al. 2002; Pescatori et al. 2007), each one batch.
- Comma-separated. Column 1 is a probe identifier, column 2 is `PROT` (gene symbol), and the remaining columns are samples whose names start with `DMD` or `NOR` (the class label).
- Probes with an empty or `!?` symbol are dropped, probes are averaged per gene, the two studies are merged on shared genes, and values are log2-transformed.

## Running a single setting

`run_BRAID.R` and `run_BEGONE.R` are a convenience for readers who want one setting at a time. They are not needed for the reproduction steps above, and they use exactly the same code paths, so a setting gives the same results either way.

**BRAID** takes four optional arguments; the defaults are shown.

```bash
Rscript run_BRAID.R VP combat standard native
#                   │  │      │        └ native | single_batch_sim
#                   │  │      └ standard (sample-wise) | batchlinked (BEAMs)
#                   │  └ combat | limma
#                   └ VP | HP
```

- `native`: the ground truth is the batch-corrected complete data.
- `single_batch_sim`: the largest native batch is split into four class-stratified pseudo-batches, known batch effects are injected (additive γ = 2.5, multiplicative δ = 0.25), and the data before injection are the ground truth.

**BEGONE** takes the dataset as its only argument.

```bash
Rscript run_BEGONE.R VP        # VP | HP
```

When sourcing a script in R instead of using `Rscript`, edit the default values at the bottom of the file. Shared settings (iteration count, K values, the BEGONE grid) live in `R/BRAID_functions.R`.

Long runs are checkpointed. If a run is interrupted, running the same command again resumes where it stopped and gives the same results as an uninterrupted run. A checkpoint written with different settings is refused; delete it to start fresh.

## Reproducing the manuscript

Nothing has to be edited in any script. Each figure script runs the experiments it needs and skips anything already finished, so the scripts can be run in any order and re-run safely.

**Main manuscript — one script:**

```bash
Rscript main_figures.R
```

This runs the eight native BRAID settings (including the parameter sweeps and the limma settings, which the Figure 3 q-values depend on) and the VP BEGONE sweep, then draws Figures 2–6.

**Supplementary — four scripts, in any order:**

```bash
Rscript supplementary/S01_batch_effect_strength.R    # simulation only; no input data needed
Rscript supplementary/S02_feature_subsampling.R      # input data only
Rscript supplementary/S03_mice_rf_limitations.R      # input data only
Rscript supplementary/S04_braid.R                    # runs all ten BRAID settings
Rscript supplementary/S05_begone.R                   # runs both BEGONE sweeps
```

In RStudio, open `BRAID.Rproj` and source any of these files; they locate the repository themselves.

BEGONE is the most expensive step: each dataset requires 100 grid cells × 10 replicates × 9 methods × 2 batch conditions = 18,000 imputations.

## Outputs

```
results/
├── braid/braid_[sim_]10iter_<dataset>_<beca>_<standardmissing|beams>/
│   ├── fast_utility.csv           one row per method × batch condition × iteration
│   ├── rmse, logfc_err, fscore    boxplots for this setting (.png and .pdf)
│   └── *_param_sweep/             parameter sweeps (VP / ComBat / sample-wise only)
├── begone/begone_10iter_<dataset>_combat_standard_9methods/
│   ├── begone_grid.csv            one row per method × grid cell × replicate
│   └── begone_heatmap_*           heatmaps for this run
├── figures/main/                  Figures 3–6
├── figures/supplementary/         Figures S1–S28
├── tables/                        statistics behind the figures, Tables S1–S2
├── batch_effect_strength/         raw values for Figures S1–S2
├── mice_rf_limitations/           raw values for Figures S4–S6
└── feature_subsampling/           raw values for Figure S3
```

**`fast_utility.csv`** columns: `method`, `set` (batch-ve or batch+ve), `RMSE`, `iteration`, `n_features_scored`, `logFC_err`, `TPR`, `FPR`, `Precision`, `Fscore`, `de_truth` (spikein or reference), and the setting (`dataset`, `beca`, `missingness`, `design`).

**`begone_grid.csv`** columns: `method`, `add_shift`, `mult_scale`, `iteration`, and for each metric its batch-ve value, batch+ve value and gap (batch+ve minus batch-ve): `batch_ve`/`batch_pve`/`gap` (RMSE), `logfc_ve`/`logfc_pve`/`logfc_gap`, `fscore_ve`/`fscore_pve`/`fscore_gap`.

Every figure is written as a PNG (350 dpi) and a vector PDF.

## Figure map

| Figure | File in `results/figures/` | Produced by | Needs |
|---|---|---|---|
| 2 | `main/Fig2_pca_*`, `main/Fig2_density_*` | `main_figures.R` | nothing |
| 3 | `main/Fig3_logratio_forest_combat` | `main_figures.R` | 8 native BRAID settings |
| 4 | `main/Fig4_begone_variance_decomposition_VP` | `main_figures.R` | BEGONE VP |
| 5A | `main/Fig5A_begone_degradation_VP` | `main_figures.R` | BEGONE VP |
| 5B | `main/Fig5B_begone_heatmap_logfc_VP` | `main_figures.R` | BEGONE VP |
| 6A | `main/Fig6A_param_sweep_KNN-sample` | `main_figures.R` | BRAID VP / ComBat / sample-wise |
| 6B | `main/Fig6B_knn_cross_batch_neighbours` | `main_figures.R` | VP input data |

Figure 2 is written as eight panels: a PCA and a density plot for each of the four combinations of additive shift (0, 2.5) and multiplicative scale (0, 0.25), each titled with its parameters.
| S1, S2 | `supplementary/FigS1*`, `FigS2*` | `S01_batch_effect_strength.R` | nothing |
| S3 | `supplementary/FigS3A*`, `FigS3B*` | `S02_feature_subsampling.R` | VP input data |
| S4–S6 | `supplementary/FigS4*` – `FigS6*` | `S03_mice_rf_limitations.R` | input data |
| S7, S8 | `supplementary/FigS7*`, `FigS8*` | `S04_braid.R` | 8 native BRAID settings |
| S9–S19 | `supplementary/FigS9*` – `FigS19*` | `S04_braid.R` | all 10 BRAID settings |
| S20 | `supplementary/FigS20A*`, `FigS20B*` | `S04_braid.R` | BRAID VP / ComBat / sample-wise |
| S21–S26 | `supplementary/FigS21*` – `FigS26*` | `S05_begone.R` | BEGONE VP and HP |
| S27, S28 | `supplementary/FigS27*`, `FigS28*` | `S05_begone.R` | BEGONE HP |
| Tables S1, S2 | `tables/TableS1*`, `tables/TableS2*` | `S05_begone.R` | BEGONE VP |

The ComBat forest plot (Figure 3) needs the limma settings too, because its q-values are Benjamini–Hochberg adjusted jointly over all settings. BRAID settings without a supplementary figure (HP with limma) are written with the prefix `Figextra`.

## Methods summary

**Missing values.** 30% of values are removed with a 3:7 MCAR:MNAR ratio (Jin et al. 2021). Features with more than 60% missing values are then dropped. In the BEAMs mode (`batchlinked`), the MNAR draw is made on batch means and applied to every sample of that batch.

**Imputation methods.**
- KNN-sample: K = 4 in BRAID for both datasets. BEGONE instead sets K from the batch structure (smallest batch-class group − 1), giving K = 6 for VP and K = 5 for HP (computed on its four pseudo-batches).
- KNN-feature: K = 10.
- SVD: rank 3.
- QRILC: `tune.sigma` = 0.5.
- NMFBatch: rank 4, 1000 iterations.
- Mean imputation.
- MICE-PMM and MICE-norm: m = 1, maxit = 5.
- RF: missForest, 100 trees.

MICE and RF use samples, not features, as variables (manuscript Section 2.6). NMFBatch corrects batch effects itself, so no further correction is applied to its output.

**Batch effect correction.** ComBat (primary) or `limma::removeBatchEffect` (secondary).

**Metrics.**
- RMSE over the imputed values.
- Median absolute logFC error: over the yeast and *E. coli* proteins for VP; over all genes for HP.
- Differential-expression F-score (limma, |logFC| threshold 0.5, FDR < 0.05). The truth is the spike-in design for VP and the differentially expressed genes called on the ground truth for HP.

**Statistical testing.** batch-ve and batch+ve are compared with paired Wilcoxon signed-rank tests on iteration-paired values, one test per method (Figures S9–S19). The log-ratio forest plots (Figures 3, S7, S8) summarise the same pairing as log((batch+ve + c) / (batch-ve + c)) per iteration, with Hodges-Lehmann estimates, exact Wilcoxon 95% confidence intervals, and Benjamini-Hochberg adjustment across every method, dataset, correction method, missingness mode and metric.

**Figure 2 illustration.** A Gaussian matrix (1000 features x 20 samples, 5 batches of 4; data seed 1234567890, simulation seed 1) shown before and after affine batch effects.

**BEGONE simulation.** Per batch, an additive shift L ~ N(α, γ/5) with α ~ U(−γ/2, γ/2), and a multiplicative scale S ~ N(β, δ/5) with β ~ U(1 − δ, 1 + δ), applied as (X + L)·S. The grid covers additive shift 0–4.5 and multiplicative scale 0–0.45, with 10 replicates per cell, each on a fresh random 35% of features.

**Reproducibility.** Every random step is seeded:
- BRAID iteration *i* uses seed 123456 + *i*.
- BEGONE cell *c*, replicate *i* uses seed 424242 + 1000·*c* + *i*.
- Pseudo-batches use seed 42.

Results therefore do not depend on which settings were run before, or whether a run was resumed from a checkpoint.

## Environment variables (optional)

| Variable | Purpose | Default |
|---|---|---|
| `BRAID_ROOT` | repository root | located automatically |
| `BRAID_DATA` | folder containing the dataset folders | `<root>/data` |
| `BRAID_RESULTS` | output folder | `<root>/results` |

## Licence

Apache License 2.0; see `LICENSE`.
