###############################################################################
## BRAID_functions.R
## Function library for BRAID (Batch-effect Reconciliation and Anomaly-aware
## Imputation for Data) and BEGONE (Batch Effect Gradient Observation on Null
## Estimation). Sourced by every script in this repository.
##
## Conventions: expression matrices are features (rows) x samples (columns),
## log2-transformed and complete before any missing values are simulated.
## Every stochastic step re-seeds itself from a deterministic seed (seed0 + i in
## BRAID, seed0 + cell * 1000 + i in BEGONE), so a setting gives identical
## results whether it is run alone or after other settings.
###############################################################################
suppressPackageStartupMessages({
  library(sva); library(limma); library(impute); library(imputeLCMD)
  library(mice); library(missForest); library(ggplot2)
})

## ---------------------------------------------------------------------------
## Paths, labels and palettes
## ---------------------------------------------------------------------------
## Input data live in <repo>/data unless BRAID_DATA points elsewhere.
data_path <- function(...) {
  base <- Sys.getenv("BRAID_DATA", unset = file.path(BRAID_ROOT, "data"))
  p <- file.path(base, ...)
  if (!file.exists(p)) stop("Input not found: ", p, "\nSee README.md, section 'Input data'.")
  p
}
results_path <- function(...) {
  p <- file.path(Sys.getenv("BRAID_RESULTS", unset = file.path(BRAID_ROOT, "results")), ...)
  dir.create(p, showWarnings = FALSE, recursive = TRUE); p
}
braid_setting_dir <- function(dataset, beca, missingness, design = "native", iterations = 10)
  results_path("braid", sprintf("braid_%s%diter_%s_%s_%s", if (design == "single_batch_sim") "sim_" else "",
                                iterations, dataset, beca, MISS_TAG[[missingness]]))
begone_run_dir <- function(dataset, iterations = 10, beca = "combat", missingness = "standard", n_methods = 9)
  results_path("begone", sprintf("begone_%diter_%s_%s_%s_%dmethods", iterations, dataset, beca, missingness, n_methods))

MVI_METHODS <- c("KNN-sample", "KNN-feature", "SVD", "QRILC", "NMFBatch", "Mean", "MICE-PMM", "MICE-norm", "RF")
MISS_TAG    <- c(standard = "standardmissing", batchlinked = "beams")
MISS_LABEL  <- c(standard = "Sample-wise", batchlinked = "BEAMs")
METRIC_LABEL <- c(RMSE = "RMSE", logFC_err = "Absolute logFC error", Fscore = "F-score")
PARAM_LABEL <- c(ksamp = "K (samples)", kfeat = "K (features)", sigma = "Sigma")
BRAID_SET_FILL   <- c("batch-ve" = "#AED9F2", "batch+ve" = "#F6B4B0", "partial" = "#F7DDA6")
BRAID_SET_COLOUR <- c("batch-ve" = "#1B6FA8", "batch+ve" = "#B23A2E", "partial" = "#B8792A")
BEGONE_PALETTE <- list(batch_pve  = c(low = "#8Fd9FB", high = "#022333", na = "grey88"),
                       logfc_pve  = c(low = "#bbe8c8", high = "#252e28", na = "grey88"),
                       fscore_pve = c(low = "#ffebb0", high = "#502F00", na = "#FFFFFF"))
BASE_SIZE <- 9
DPI <- 350

thin_theme <- function(legend = "bottom")
  theme_bw(base_size = BASE_SIZE) +
  theme(legend.position = legend, legend.key.size = unit(0.8, "lines"),
        legend.box.spacing = unit(4, "pt"), panel.grid.minor = element_blank(),
        plot.title = element_text(size = BASE_SIZE + 0.5, face = "bold"),
        plot.subtitle = element_text(size = BASE_SIZE - 1, colour = "grey30"),
        plot.margin = margin(4, 5, 2, 3))

## Write a figure as PNG (DPI) and vector PDF.
save_fig <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = DPI, limitsize = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE,
         device = if (capabilities("cairo")) cairo_pdf else "pdf")   # cairo keeps non-Latin-1 glyphs such as the em dash
  invisible(p)
}

## ---------------------------------------------------------------------------
## Dataset loaders. Each returns list(df, batch_f, class_f), complete and log2.
## ---------------------------------------------------------------------------
## Van Puyvelde: DDA runs from four instruments (one batch each), proteins
## present in every run.
load_VP <- function() {
  dir <- data_path("Van Puyvelde dataset")
  inst <- c("HYE5600735", "HYE6600735", "HYEqe735", "HYEtims735")
  design <- lapply(inst, function(f) read.table(file.path(dir, paste0(f, "_LFQ_FragPipe_design.tsv")), sep = "\t", header = TRUE))
  ints   <- lapply(inst, function(f) read.table(file.path(dir, paste0(f, "_LFQ_FragPipe_pro_intensity.tsv")), sep = "\t", header = TRUE))
  ints   <- lapply(ints, function(d) { m <- d[, c(1, 3:ncol(d))]; rownames(m) <- d$Protein; m })
  comb <- Reduce(function(x, y) merge(x, y, by = "Protein", all = TRUE), ints)
  comb <- na.omit(comb); comb[comb == 0] <- NA; comb <- na.omit(comb)
  rownames(comb) <- comb$Protein
  list(df = log2(comb[, -1]),
       batch_f = as.factor(rep(1:4, times = sapply(design, nrow))),
       class_f = as.factor(as.numeric(as.factor(do.call(rbind, design)$condition))))
}

## Haslett & Pescatori: two DMD microarray studies (one batch each), probes
## collapsed to gene symbols by averaging, merged on shared genes.
load_HP <- function() {
  dir <- data_path("Haslett & Pescatori datasets")
  a <- read.csv(file.path(dir, "DMD-HaslettData.csv")); b <- read.csv(file.path(dir, "DMD-PescatoriData.csv"))
  ca <- substr(colnames(a)[3:ncol(a)], 1, 3); cb <- substr(colnames(b)[3:ncol(b)], 1, 3)   # DMD / NOR
  a[a$PROT %in% c("", "!?"), "PROT"] <- NA; a <- na.omit(a)
  b[b$PROT %in% c("", "!?"), "PROT"] <- NA; b <- na.omit(b)
  sc <- function(d) grep("^(DMD|NOR)", colnames(d), value = TRUE)
  s1 <- aggregate(a[, sc(a)], by = list(PROT = a$PROT), FUN = mean, na.rm = TRUE)
  s2 <- aggregate(b[, sc(b)], by = list(PROT = b$PROT), FUN = mean, na.rm = TRUE)
  cm <- merge(s1, s2, by = "PROT"); rownames(cm) <- cm$PROT
  ex <- as.matrix(log2(cm[, -1])); ex[!is.finite(ex)] <- NA; ex <- na.omit(ex)
  list(df = as.data.frame(ex),
       batch_f = as.factor(c(rep(1, ncol(a[, -c(1, 2)])), rep(2, ncol(b[, -c(1, 2)])))),
       class_f = as.factor(c(ca, cb)))
}
load_dataset <- function(name) switch(name, VP = load_VP(), HP = load_HP(), stop("Unknown dataset: ", name))

## Class-stratified pseudo-batch labels (every pseudo-batch keeps every class).
assign_pseudo_batches <- function(class_f, n_batches, seed) {
  set.seed(seed)
  pseudo <- integer(length(class_f))
  for (cl in levels(class_f)) {
    idx <- which(class_f == cl)
    pseudo[idx] <- sample(rep(seq_len(n_batches), length.out = length(idx)))
  }
  as.factor(pseudo)
}

## Single-batch validation design: take the largest native batch (genuinely
## batch-free), split it into pseudo-batches and inject known affine batch
## effects. The pre-injection matrix is returned as the ground truth.
make_single_batch_sim <- function(d, n_pseudo = 4, shift = 2.5, scale = 0.25, seed = 42) {
  chosen <- names(sort(table(d$batch_f), decreasing = TRUE))[1]
  keep <- as.character(d$batch_f) == chosen
  if (sum(keep) < n_pseudo * 2) stop("Batch ", chosen, " is too small to split into ", n_pseudo, " pseudo-batches.")
  single <- d$df[, keep, drop = FALSE]; cls <- droplevels(d$class_f[keep])
  pseudo <- assign_pseudo_batches(cls, n_pseudo, seed)
  sim <- sim_batch_effects(as.matrix(single), pseudo, shift = shift, scale = scale, seed = seed)
  cat(sprintf("  single-batch design: batch %s (%d samples) -> %d pseudo-batches, shift = %.2f, scale = %.2f\n",
              chosen, sum(keep), n_pseudo, shift, scale))
  list(df = as.data.frame(sim), batch_f = pseudo, class_f = cls, external_truth = single)
}

## ---------------------------------------------------------------------------
## Missing-value simulation (Jin et al. 2021: MCAR + MNAR mixture)
## ---------------------------------------------------------------------------
## mv_prop = c(total MV fraction, MNAR fraction). With simbeam = TRUE the MNAR
## draw is made on batch-level means and broadcast to every sample of that
## batch (batch-effect-associated missingness, BEAMs).
simulate_mvs <- function(data, mv_prop, simbeam = FALSE, batch_factor = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Actual_data <- data
  if (simbeam) {
    batches <- unique(batch_factor); beams_data <- c()
    for (i in batches) beams_data <- cbind(beams_data, apply(data[, which(batch_factor == i)], 1, mean, na.rm = TRUE))
    colnames(beams_data) <- seq(length(batches)); data <- beams_data
  }
  D <- as.matrix(data); n <- length(D); alpha <- mv_prop[1]; beta <- mv_prop[2]
  tmat <- matrix(NA, ncol = ncol(D), nrow = nrow(D))
  for (col in 1:ncol(tmat)) tmat[, col] <- rnorm(nrow(tmat), mean = quantile(data[, col], alpha, na.rm = TRUE), sd = 0.3)
  pmat <- matrix(rbinom(n, 1, beta), ncol = ncol(D), nrow = nrow(D))
  indicator <- (D < tmat) * 1 + pmat
  msdata <- D; msdata[indicator == 2] <- NA
  if (simbeam) {
    msdata <- Actual_data
    for (bi in seq_along(batches)) msdata[which(indicator[, bi] == 2), which(batch_factor == batches[bi])] <- NA
  }
  total_mcar <- n * alpha * (1 - beta)
  mcar_inds <- sample(which(is.finite(msdata)), round(total_mcar), replace = FALSE)
  if (simbeam) {
    mcar_inds <- sample(which(is.finite(indicator)), round(total_mcar), replace = FALSE)
    indicator[mcar_inds] <- 3
    for (bi in seq_along(batches)) msdata[which(indicator[, bi] == 3), which(batch_factor == batches[bi])] <- NA
  } else msdata[mcar_inds] <- NA
  list(mcar_ind = mcar_inds, mnar_ind = indicator, msdata = msdata)
}
simulate_mvs_mode <- function(data, mv_prop, batch_factor, missingness = c("standard", "batchlinked"), seed = NULL)
  simulate_mvs(data, mv_prop, simbeam = match.arg(missingness) == "batchlinked", batch_factor = batch_factor, seed = seed)

## Drop features whose missingness exceeds `threshold` (BRAID Step 3).
filter_high_mv_features <- function(data, threshold = 0.6)
  data[apply(data, 1, function(x) mean(is.na(x)) <= threshold), , drop = FALSE]

## Copy an MV mask onto the same features of another matrix (BRAID Step 5).
transfer_mv_pattern <- function(full, mask) { s <- full[rownames(mask), , drop = FALSE]; s[is.na(mask)] <- NA; s }

## Random subset of a fraction of features (rows).
subsample_features <- function(mat, frac, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mat[sort(sample.int(nrow(mat), size = ceiling(nrow(mat) * frac))), , drop = FALSE]
}

## ---------------------------------------------------------------------------
## Batch-effect simulation (BEGONE Step 2)
## ---------------------------------------------------------------------------
## Per batch i: alpha_i ~ U(-shift/2, shift/2), L_ig ~ N(alpha_i, shift/5);
## beta_i ~ U(1 - scale, 1 + scale), S_ig ~ N(beta_i, scale/5).
## Applied as X_batch = (X + L) * S. A zero shift (scale) skips that component.
sim_batch_effects <- function(data, batch_assignments, shift = 0.5, scale = 0.1, seed = 1234567890) {
  if (!is.null(seed)) set.seed(seed)
  data <- as.matrix(data); n_feat <- nrow(data)
  for (g in unique(batch_assignments)) {
    idx <- which(batch_assignments == g); L <- rep(0, n_feat); S <- rep(1, n_feat)
    if (shift != 0) { alpha <- runif(1, -(shift / 2), shift / 2); L <- rnorm(n_feat, mean = alpha, sd = shift / 5) }
    if (scale != 0) { beta <- runif(1, 1 - scale, 1 + scale); S <- rnorm(n_feat, mean = beta, sd = scale / 5) }
    data[, idx] <- (data[, idx, drop = FALSE] + L) * S
  }
  data
}

## ---------------------------------------------------------------------------
## Batch-effect correction (ComBat primary, limma::removeBatchEffect secondary)
## ---------------------------------------------------------------------------
apply_beca <- function(data, batch, method = c("combat", "limma"), class = NULL) {
  method <- match.arg(method); batch <- as.factor(batch); dat <- as.matrix(data)
  if (method == "combat") return(sva::ComBat(dat = .stabilize_for_combat(dat, batch), batch = batch))
  limma::removeBatchEffect(dat, batch = batch, design = if (is.null(class)) NULL else model.matrix(~ as.factor(class)))
}
## Constant-fill imputers can leave a feature with zero within-batch variance,
## which breaks ComBat. Such (feature x batch) blocks receive a deterministic,
## negligible (1e-4 x feature SD) mean-zero ramp; all other blocks are untouched.
.stabilize_for_combat <- function(data, batch, eps_frac = 1e-4) {
  gsd <- stats::sd(as.vector(data), na.rm = TRUE); if (!is.finite(gsd) || gsd == 0) gsd <- 1
  for (b in levels(batch)) {
    cc <- which(batch == b); if (length(cc) < 2) next
    rv <- apply(data[, cc, drop = FALSE], 1, stats::var, na.rm = TRUE)
    bad <- which(!is.finite(rv) | rv <= .Machine$double.eps)
    if (length(bad)) {
      fsd <- apply(data[bad, , drop = FALSE], 1, stats::sd, na.rm = TRUE); fsd[!is.finite(fsd) | fsd == 0] <- gsd
      data[bad, cc] <- data[bad, cc] + outer(eps_frac * fsd, as.numeric(scale(seq_along(cc))))
    }
  }
  data
}
create_ground_truth   <- function(data, batch_factor, beca = "combat", class = NULL) apply_beca(data, batch_factor, beca, class)
correct_batch_effects <- function(imputed_list, batch_factor, beca = "combat", class = NULL)
  lapply(imputed_list, apply_beca, batch = batch_factor, method = beca, class = class)

## ---------------------------------------------------------------------------
## Imputation
## ---------------------------------------------------------------------------
## KNN-sample: each missing value is the inverse-distance weighted mean of the
## K nearest samples (euclidean, over observed features) observed for that feature.
do.knn <- function(df, k) {
  df <- as.matrix(df); imp <- df
  dist_mat <- as.matrix(dist(t(df), method = "euclidean"))
  for (i in unique(which(is.na(df), arr.ind = TRUE)[, 1])) {
    obs <- which(!is.na(df[i, ]))
    for (s in which(is.na(df[i, ]))) {
      d <- dist_mat[s, obs]
      if (length(d) < k) stop("Fewer than K finite distances found")
      nn <- order(d)[1:k]; w <- 1 / d[nn] / sum(1 / d[nn])
      imp[i, s] <- w %*% df[i, obs[nn]]
    }
  }
  imp
}

## Samples are the variables for MICE and RF (see manuscript Section 2.6);
## samples with no observed value are excluded from fitting and returned as is.
.fit_on_samples <- function(mat, fun) {
  keep <- colSums(!is.na(mat)) > 0
  if (any(keep)) mat[, keep] <- as.matrix(fun(mat[, keep, drop = FALSE]))
  mat
}
## A method that errors leaves its matrix unimputed and raises a warning.
impute_one <- function(mat, method, ksamp = 4, kfeat = 10, ntrees = 100, sigma = 0.5) {
  mat <- as.matrix(mat)
  res <- tryCatch(switch(method,
    "MICE-PMM"    = .fit_on_samples(mat, function(x) mice::complete(mice::mice(as.data.frame(x), m = 1, maxit = 5, printFlag = FALSE))),
    "MICE-norm"   = .fit_on_samples(mat, function(x) mice::complete(mice::mice(as.data.frame(x), meth = "norm", m = 1, maxit = 5, printFlag = FALSE))),
    "RF"          = .fit_on_samples(mat, function(x) missForest(x, ntree = ntrees)$ximp),
    "KNN-feature" = impute.knn(mat, k = min(kfeat, max(1L, nrow(mat) - 1L)), rowmax = 0.99, colmax = 0.99)$data,
    "KNN-sample"  = do.knn(mat, k = min(ksamp, max(1L, ncol(mat) - 1L))),
    "SVD"         = impute.wrapper.SVD(mat, K = min(3, min(dim(mat)))),
    "QRILC"       = impute.QRILC(mat, tune.sigma = sigma)[[1]],
    "Mean"        = { v <- rowMeans(mat, na.rm = TRUE); v[is.nan(v)] <- NA; i <- which(is.na(mat)); mat[i] <- v[row(mat)[i]]; mat },
    stop("Unknown imputation method: ", method)),
    error = function(e) { warning("Could not apply ", method, "; leaving it unimputed. Reason: ", conditionMessage(e)); mat })
  res <- as.matrix(res); dim(res) <- dim(mat); dimnames(res) <- dimnames(mat); res
}
## Methods run in the order given; the order matters for reproducibility
## because MICE, RF and QRILC draw from the shared random-number stream.
impute_methods <- function(data, methods, verbose = TRUE, ...) {
  out <- list()
  for (m in methods) {
    tt <- system.time(out[[m]] <- impute_one(data, m, ...))["elapsed"]
    if (verbose) message(sprintf("Applying method: %-12s | Time taken: %.2f seconds", m, tt))
  }
  out
}

## NMFBatch: one-step imputation + batch correction (rank k, GAM batch removal
## with full model batch + class vs reduced model class). Its output is already
## batch-corrected, so no BECA is applied on top. Its data preparation can drop
## features, so all scoring re-aligns on shared rownames.
run_nmfbatch <- function(mat, batch_factor, class_factor, k = 4, max_iter = 1000) {
  suppressPackageStartupMessages(library(NMFBatch))
  mat <- as.matrix(mat)
  meta <- data.frame(batch = as.factor(batch_factor), class = as.factor(class_factor), row.names = colnames(mat))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(Uncorrected = mat), colData = meta)
  se_prep <- prepare_data(se, k = k)
  fit <- batch_factorization(se_prep, k = k, max.iter = max_iter)
  fit_c <- gam_batch(NMF_fit = fit, se = se_prep, k = k, full_formula = "batch + class", reduced_formula = "class")
  out <- reconstruct_with_residuals(se = se_prep, original_NMF_fit = fit, corrected_NMF_fit = fit_c,
                                    impute = TRUE)@assays@data@listData$NMFcorrected
  dropped <- setdiff(rownames(mat), rownames(out))
  if (length(dropped)) message(sprintf("  [NMFBatch] dropped %d/%d input features", length(dropped), nrow(mat)))
  out
}
## Runs NMFBatch on the batch-ve and batch+ve matrices (in that order) and
## splices the results into the two method lists.
.add_nmfbatch <- function(ve_list, pve_list, mis_ve, mis_pve, batch_factor, class_factor, k, max_iter) {
  try_nmf <- function(m, lab) tryCatch(run_nmfbatch(m, batch_factor, class_factor, k, max_iter),
                                       error = function(e) { message("  [NMFBatch/", lab, "] failed: ", conditionMessage(e)); NULL })
  a <- try_nmf(mis_ve, "batch-ve"); b <- try_nmf(mis_pve, "batch+ve")
  if (!is.null(a)) ve_list[["NMFBatch"]] <- a
  if (!is.null(b)) pve_list[["NMFBatch"]] <- b
  list(ve = ve_list, pve = pve_list)
}

## Keep only the requested methods whose packages are installed.
available_methods <- function(methods) {
  pkg <- c("NMFBatch" = "NMFBatch")
  for (m in intersect(methods, names(pkg)))
    if (!requireNamespace(pkg[[m]], quietly = TRUE)) { cat("  note: dropping", m, "(package", pkg[[m]], "not installed)\n"); methods <- setdiff(methods, m) }
  methods
}

## KNN-sample K: one less than the smallest batch-class group (minimum 1).
batch_aware_k_sample <- function(batch_factor, class_factor) max(1L, min(table(as.factor(batch_factor), as.factor(class_factor))) - 1L)

## ---------------------------------------------------------------------------
## Evaluation metrics
## ---------------------------------------------------------------------------
## RMSE over the cells that are missing in `mis`.
Rmse <- function(imp, mis, true) {
  idx <- which(is.na(as.matrix(mis)))
  sqrt(mean((as.matrix(imp)[idx] - as.matrix(true)[idx])^2))
}

## limma fit with FDR-adjusted q-values; applyfc = TRUE tests against a
## log fold-change threshold (limma::treat).
limma_fit <- function(df, class, applyfc = TRUE, thres = 0.5) {
  design <- model.matrix(~ as.factor(class)); colnames(design) <- c("Intercept", "Contrast")
  fit <- eBayes(lmFit(df, design = design))
  if (applyfc) fit <- treat(fit, lfc = thres)
  fit$q.value <- p.adjust(fit$p.value[, "Contrast"], method = "fdr"); fit
}
## output = "significant": features with q < 0.05; "logfc": estimated logFC;
## "DEA": TPR/FPR/Precision/Fscore against spike-in truth (YEAST/ECOLI are
## differential, HUMAN is not).
limmaDEA <- function(df, class_factor, applyfc = TRUE, set.thres = 0.5, output = c("significant", "logfc", "DEA")) {
  output <- match.arg(output)
  fit <- limma_fit(df, class_factor, applyfc, set.thres)
  if (output == "logfc") return(fit$coefficients[, 2])
  idx <- which(fit$q.value < 0.05); pos <- rownames(df)[idx]
  if (output == "significant") return(pos)
  neg <- if (length(idx)) rownames(df)[-idx] else character(0)
  TP <- sum(grepl("YEAST", pos)) + sum(grepl("ECOLI", pos)); FP <- sum(grepl("HUMAN", pos))
  TN <- sum(grepl("HUMAN", neg)); FN <- sum(grepl("YEAST", neg)) + sum(grepl("ECOLI", neg))
  TPR <- TP / (TP + FN); Precision <- TP / (TP + FP)
  c(TPR = TPR, FPR = FP / (FP + TN), Precision = Precision, Fscore = 2 * Precision * TPR / (Precision + TPR))
}

## Median absolute logFC error against the truth, over the spike-in
## differential features (all features when no spike-in tags exist, e.g. HP).
logfc_error_one <- function(imputed, true, class_factor, dap_pattern = "YEAST|ECOLI") {
  dap <- rownames(true)[grepl(dap_pattern, rownames(true), ignore.case = TRUE)]
  if (!length(dap)) dap <- rownames(true)
  getlfc <- function(x) tryCatch({ v <- limmaDEA(x, class_factor, applyfc = FALSE, output = "logfc")
                                   setNames(as.numeric(v), rownames(x)[seq_along(v)]) }, error = function(e) NULL)
  tl <- getlfc(true); il <- getlfc(imputed)
  if (is.null(tl) || is.null(il)) return(NA_real_)
  dap <- intersect(dap, intersect(names(tl), names(il)))
  if (!length(dap)) return(NA_real_)
  stats::median(abs(il[dap] - tl[dap]), na.rm = TRUE)
}

## Differential-expression recovery. "spikein": truth from YEAST/ECOLI/HUMAN
## tags (VP). "reference": truth = DE set called on the ground-truth matrix (HP).
evaluate_downstream <- function(imputed, class_factor, truth_source = c("spikein", "reference"),
                                reference_sig = NULL, universe = NULL, applyfc = TRUE, set.thres = 0.5) {
  truth_source <- match.arg(truth_source)
  if (truth_source == "spikein") {
    m <- limmaDEA(as.matrix(imputed), class_factor, applyfc, set.thres, output = "DEA")
    return(data.frame(TPR = unname(m["TPR"]), FPR = unname(m["FPR"]), Precision = unname(m["Precision"]), Fscore = unname(m["Fscore"])))
  }
  pred <- tryCatch(limmaDEA(as.matrix(imputed), class_factor, applyfc, set.thres, output = "significant"), error = function(e) character(0))
  universe <- unique(universe); ref <- intersect(reference_sig, universe); pred <- intersect(pred, universe)
  TP <- length(intersect(pred, ref)); FP <- length(setdiff(pred, ref)); FN <- length(setdiff(ref, pred))
  TN <- length(universe) - length(union(pred, ref))
  Precision <- if (TP + FP > 0) TP / (TP + FP) else NA_real_; TPR <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  Fscore <- if (!is.na(Precision) && !is.na(TPR) && Precision + TPR > 0) 2 * Precision * TPR / (Precision + TPR) else NA_real_
  data.frame(TPR = TPR, FPR = if (FP + TN > 0) FP / (FP + TN) else NA_real_, Precision = Precision, Fscore = Fscore)
}
.de_truth_source <- function(rn) if (any(grepl("YEAST|ECOLI|HUMAN", rn, ignore.case = TRUE))) "spikein" else "reference"
.ref_sig <- function(truth, class_factor, applyfc = TRUE, thres = 0.5)
  tryCatch(limmaDEA(truth, class_factor, applyfc, thres, output = "significant"), error = function(e) character(0))

## ---------------------------------------------------------------------------
## Checkpointing: a long run resumes from its last completed unit. A checkpoint
## written with different settings is refused (delete it to start fresh).
## ---------------------------------------------------------------------------
.ckpt_read <- function(path, meta) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ck <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(ck) && !identical(ck$meta, meta)) stop("Checkpoint '", path, "' was written with different settings; delete it to start fresh.")
  ck
}
.ckpt_write <- function(path, obj) if (!is.null(path)) { saveRDS(obj, paste0(path, ".tmp")); file.rename(paste0(path, ".tmp"), path) }

## ---------------------------------------------------------------------------
## BRAID
## ---------------------------------------------------------------------------
## Ground truth (Step 2) = BECA(complete data), or `external_truth` in the
## single-batch design. Per iteration i (seed seed0 + i):
##   batch-ve = impute(ground truth + MVs)                          (Steps 3, 6)
##   batch+ve = BECA(impute(raw data + the same MVs))                (Steps 4-7)
## and both are scored against the ground truth (Step 8): RMSE, median
## absolute logFC error, and DE-recovery F-score.
run_braid <- function(df, batch_f, class_f, mv_prop, iterations, mvi_methods = MVI_METHODS,
                      missingness = "standard", beca = "combat", ksamp = 4, kfeat = 10, ntrees = 100, sigma = 0.5,
                      external_truth = NULL, applyfc = TRUE, set.thres = 0.5, mv_filter = 0.6,
                      nmf_k = 4, nmf_max_iter = 1000, seed0 = 123456, checkpoint = NULL, verbose = TRUE) {
  run_nmf <- "NMFBatch" %in% mvi_methods; mvi_methods <- setdiff(mvi_methods, "NMFBatch")
  raw <- as.matrix(df); batch_factor <- as.factor(batch_f); class_factor <- as.factor(class_f)
  stopifnot(length(batch_factor) == ncol(raw), length(class_factor) == ncol(raw))
  truth_source <- .de_truth_source(rownames(raw))
  GT <- if (is.null(external_truth)) create_ground_truth(raw, batch_factor, beca, class_factor) else as.matrix(external_truth)
  stopifnot(identical(dim(GT), dim(raw)), identical(rownames(GT), rownames(raw)))
  ref_sig_full <- if (truth_source == "reference") .ref_sig(GT, class_factor, applyfc, set.thres) else NULL
  imp_args <- list(ksamp = ksamp, kfeat = kfeat, ntrees = ntrees, sigma = sigma)
  meta <- list(iterations = iterations, seed0 = seed0, methods = mvi_methods, run_nmf = run_nmf, beca = beca,
               missingness = missingness, dim = dim(raw), external = !is.null(external_truth), imp_args = imp_args)
  ck <- .ckpt_read(checkpoint, meta)
  save <- if (is.null(ck)) NULL else ck$save; start_i <- if (is.null(ck)) 1L else ck$next_i
  for (i in seq_len(iterations)[seq_len(iterations) >= start_i]) {
    if (verbose) cat("  iteration", i, "of", iterations, "\n")
    mv    <- simulate_mvs_mode(GT, mv_prop, batch_factor, missingness, seed = seed0 + i)
    misA  <- filter_high_mv_features(mv$msdata, mv_filter)
    truth <- GT[rownames(misA), , drop = FALSE]
    misB  <- transfer_mv_pattern(raw, misA)
    A_imp <- do.call(impute_methods, c(list(misA, mvi_methods), imp_args))
    B_imp <- do.call(impute_methods, c(list(misB, mvi_methods), imp_args))
    B_cor <- correct_batch_effects(B_imp, batch_factor, beca, class_factor)
    if (run_nmf) { nm <- .add_nmfbatch(A_imp, B_cor, misA, misB, batch_factor, class_factor, nmf_k, nmf_max_iter); A_imp <- nm$ve; B_cor <- nm$pve }
    for (m in names(A_imp)) for (pth in c("batch-ve", "batch+ve")) {
      mat <- if (pth == "batch-ve") A_imp[[m]] else B_cor[[m]]
      if (is.null(mat)) next
      mis <- if (pth == "batch-ve") misA else misB
      rn <- intersect(rownames(mat), rownames(truth))
      mat <- mat[rn, , drop = FALSE]; tr <- truth[rn, , drop = FALSE]
      row <- data.frame(method = m, set = pth, RMSE = Rmse(mat, mis[rn, , drop = FALSE], tr), iteration = i,
                        n_features_scored = length(rn),
                        logFC_err = tryCatch(logfc_error_one(mat, tr, class_factor), error = function(e) NA_real_))
      du <- tryCatch(evaluate_downstream(mat, class_factor, truth_source, intersect(ref_sig_full, rn), rn, applyfc, set.thres),
                     error = function(e) data.frame(TPR = NA_real_, FPR = NA_real_, Precision = NA_real_, Fscore = NA_real_))
      save <- rbind(save, cbind(row, du, de_truth = truth_source))
    }
    .ckpt_write(checkpoint, list(save = save, next_i = i + 1L, meta = meta))
  }
  save
}

## Runs BRAID once per value of one method's hyperparameter.
braid_param_sweep <- function(df, batch_f, class_f, mv_prop, iterations, method, param_name, param_values,
                              checkpoint_dir = NULL, ...) {
  do.call(rbind, lapply(param_values, function(v) {
    cat(sprintf("== %s: %s = %s ==\n", method, param_name, v))
    args <- c(list(df = df, batch_f = batch_f, class_f = class_f, mv_prop = mv_prop, iterations = iterations,
                   mvi_methods = method), list(...))
    args[[param_name]] <- v
    if (!is.null(checkpoint_dir)) args$checkpoint <- file.path(checkpoint_dir, sprintf("sweep_%s_%s_%s.rds", method, param_name, v))
    res <- do.call(run_braid, args); res$param_name <- param_name; res$param_value <- v; res
  }))
}

## Collects every results/braid/<setting>/fast_utility.csv into one table.
load_braid_results <- function() {
  f <- list.files(results_path("braid"), pattern = "^fast_utility\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("No BRAID results under ", results_path("braid"), ". Run run_BRAID.R first (see README).")
  do.call(rbind, lapply(f, read.csv, stringsAsFactors = FALSE))
}

## ---------------------------------------------------------------------------
## BEGONE
## ---------------------------------------------------------------------------
## For every (additive shift, multiplicative scale) cell and iteration
## (seed seed0 + cell * 1000 + i): subsample features, inject batch effects
## into the batch-free reference, simulate MVs, then compare
## impute(reference + MVs) with BECA(impute(batch-affected + same MVs)).
run_begone <- function(clean_df, batch_f, class_f, add_levels, mult_levels, mv_prop, iterations,
                       mvi_methods = MVI_METHODS, missingness = "standard", beca = "combat",
                       ksamp = 2, kfeat = 10, ntrees = 100, sigma = 0.5, feature_subsample = FALSE,
                       fscore_applyfc = TRUE, fscore_thres = 0.5, mv_filter = 0.6,
                       nmf_k = 4, nmf_max_iter = 1000, seed0 = 424242, checkpoint = NULL, verbose = TRUE) {
  run_nmf <- "NMFBatch" %in% mvi_methods; mvi_methods <- setdiff(mvi_methods, "NMFBatch")
  G <- as.matrix(clean_df); batch_factor <- as.factor(batch_f); class_factor <- as.factor(class_f)
  truth_source <- .de_truth_source(rownames(G))
  imp_args <- list(ksamp = ksamp, kfeat = kfeat, ntrees = ntrees, sigma = sigma)
  grid <- expand.grid(add_shift = add_levels, mult_scale = mult_levels)
  meta <- list(iterations = iterations, seed0 = seed0, methods = mvi_methods, run_nmf = run_nmf, beca = beca,
               missingness = missingness, add = add_levels, mult = mult_levels, fs = feature_subsample,
               dim = dim(G), imp_args = imp_args)
  ck <- .ckpt_read(checkpoint, meta)
  out <- if (is.null(ck)) NULL else ck$out; start_gi <- if (is.null(ck)) 1L else ck$next_gi
  score <- function(mat, mis, truth, sig) {
    rn <- intersect(rownames(mat), rownames(truth))
    if (!length(rn)) return(c(NA_real_, NA_real_, NA_real_, 0))
    m <- mat[rn, , drop = FALSE]; tr <- truth[rn, , drop = FALSE]
    f <- tryCatch(evaluate_downstream(m, class_factor, truth_source, intersect(sig, rn), rn, fscore_applyfc, fscore_thres)$Fscore,
                  error = function(e) NA_real_)
    c(Rmse(m, mis[rn, , drop = FALSE], tr), logfc_error_one(m, tr, class_factor), f, length(rn))
  }
  for (gi in seq_len(nrow(grid))[seq_len(nrow(grid)) >= start_gi]) {
    a <- grid$add_shift[gi]; s <- grid$mult_scale[gi]
    if (verbose) cat(sprintf("  cell %d/%d  add = %.2f  mult = %.2f\n", gi, nrow(grid), a, s))
    for (i in seq_len(iterations)) {
      seed_i <- seed0 + gi * 1000L + i
      G_i  <- if (!isFALSE(feature_subsample)) subsample_features(G, feature_subsample, seed = seed_i) else G
      B    <- sim_batch_effects(G_i, batch_factor, shift = a, scale = s, seed = seed_i)
      mv   <- simulate_mvs_mode(G_i, mv_prop, batch_factor, missingness, seed = seed_i)
      misG <- filter_high_mv_features(mv$msdata, mv_filter)
      truth <- G_i[rownames(misG), , drop = FALSE]
      misB <- B[rownames(misG), , drop = FALSE]; misB[is.na(misG)] <- NA
      sig  <- if (truth_source == "reference") .ref_sig(truth, class_factor, fscore_applyfc, fscore_thres) else NULL
      Av <- do.call(impute_methods, c(list(misG, mvi_methods), imp_args))
      Bi <- do.call(impute_methods, c(list(misB, mvi_methods), imp_args))
      Bc <- correct_batch_effects(Bi, batch_factor, beca, class_factor)
      if (run_nmf) { nm <- .add_nmfbatch(Av, Bc, misG, misB, batch_factor, class_factor, nmf_k, nmf_max_iter); Av <- nm$ve; Bc <- nm$pve }
      for (m in names(Av)) {
        v <- score(Av[[m]], misG, truth, sig)
        p <- if (m %in% names(Bc)) score(Bc[[m]], misB, truth, sig) else c(NA_real_, NA_real_, NA_real_, 0)
        out <- rbind(out, data.frame(method = m, add_shift = a, mult_scale = s, batch_ve = v[1], batch_pve = p[1],
                                     gap = p[1] - v[1], iteration = i, n_features_ve = v[4], n_features_pve = p[4],
                                     logfc_ve = v[2], logfc_pve = p[2], logfc_gap = p[2] - v[2],
                                     fscore_ve = v[3], fscore_pve = p[3], fscore_gap = p[3] - v[3]))
      }
    }
    .ckpt_write(checkpoint, list(out = out, next_gi = gi + 1L, meta = meta))
  }
  out
}

## BEGONE statistics (manuscript Section 2.7.4), with simulation replicate
## (`iteration`) as the clustering unit.
## (1) Type II variance decomposition of outcome ~ method + add_shift + mult_scale + replicate.
begone_variance_decomposition <- function(df, outcome) {
  sub <- df[!is.na(df[[outcome]]), , drop = FALSE]
  sub$iteration <- as.factor(sub$iteration); sub$method <- as.factor(sub$method)
  aov <- as.data.frame(car::Anova(lm(as.formula(paste(outcome, "~ method + add_shift + mult_scale + iteration")), data = sub), type = 2))
  aov$term <- rownames(aov); aov$pct_var <- 100 * aov[["Sum Sq"]] / sum(aov[["Sum Sq"]]); aov$outcome <- outcome
  aov <- aov[order(-aov$pct_var), c("term", "Sum Sq", "Df", "F value", "Pr(>F)", "pct_var", "outcome")]
  rownames(aov) <- NULL; aov
}
## (2) outcome ~ method * (add_shift + mult_scale): per-method slopes with 95%
## cluster-robust (by replicate) confidence intervals.
begone_slope_model <- function(df, outcome, conf_level = 0.95) {
  sub <- df[!is.na(df[[outcome]]), , drop = FALSE]
  sub$iteration <- as.factor(sub$iteration); sub$method <- as.factor(sub$method)
  ref <- levels(sub$method)[1]
  m <- lm(as.formula(paste(outcome, "~ method * (add_shift + mult_scale)")), data = sub)
  V <- sandwich::vcovCL(m, cluster = sub$iteration); b <- coef(m); z <- qnorm(1 - (1 - conf_level) / 2)
  est_se <- function(meth, base) {   # estimate and robust SE of a method's intercept or slope
    v <- setNames(rep(0, length(b)), names(b)); v[base] <- 1
    if (meth != ref) { nm <- paste0("method", meth, if (base == "(Intercept)") "" else paste0(":", base)); if (nm %in% names(v)) v[nm] <- 1 }
    nm <- intersect(names(v), rownames(V))
    c(sum(b[v == 1]), sqrt(max(as.numeric(t(v[nm]) %*% V[nm, nm] %*% v[nm]), 0)))
  }
  tab <- do.call(rbind, lapply(levels(sub$method), function(meth) {
    i0 <- est_se(meth, "(Intercept)"); sa <- est_se(meth, "add_shift"); sm <- est_se(meth, "mult_scale")
    data.frame(method = meth, intercept = i0[1],
               slope_add_shift = sa[1], slope_add_shift_lo = sa[1] - z * sa[2], slope_add_shift_hi = sa[1] + z * sa[2],
               slope_mult_scale = sm[1], slope_mult_scale_lo = sm[1] - z * sm[2], slope_mult_scale_hi = sm[1] + z * sm[2],
               pred_min_severity = i0[1],
               pred_max_severity = i0[1] + sa[1] * max(sub$add_shift) + sm[1] * max(sub$mult_scale))
  }))
  tab <- tab[order(tab$slope_mult_scale), ]; rownames(tab) <- NULL
  list(model = m, table = tab, outcome = outcome)
}
BEGONE_OUTCOMES <- c(batch_pve = "RMSE", logfc_pve = "logFC error", fscore_pve = "Fscore")
begone_quant <- function(grid)
  list(variance = setNames(lapply(names(BEGONE_OUTCOMES), begone_variance_decomposition, df = grid), names(BEGONE_OUTCOMES)),
       slopes   = setNames(lapply(names(BEGONE_OUTCOMES), begone_slope_model, df = grid), names(BEGONE_OUTCOMES)))

plot_variance_decomposition <- function(var_list) {
  df <- do.call(rbind, var_list); df$term[df$term == "Residuals"] <- "Residual"
  lv <- c("method", "add_shift", "mult_scale", "iteration", "Residual")
  df$term <- factor(df$term, levels = lv[lv %in% df$term])
  ggplot(df, aes(outcome, pct_var, fill = term)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(pct_var >= 4, sprintf("%.0f%%", pct_var), "")), position = position_stack(vjust = 0.5), size = 3, colour = "white") +
    coord_flip() + scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "% of variance explained", fill = "Term", title = "Variance decomposition of BEGONE grid outcomes") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold"), legend.position = "bottom")
}
## Full-range degradation: slope x axis range; F-score sign-reversed so that
## positive is always worse. Methods ordered by F-score degradation (multiplicative).
plot_begone_degradation <- function(slopes, add_range, mult_range) {
  eff <- do.call(rbind, lapply(names(slopes), function(o) {
    t <- slopes[[o]]$table; met <- BEGONE_OUTCOMES[[o]]; sg <- if (met == "Fscore") -1 else 1
    rbind(data.frame(method = t$method, metric = met, axis = "mult", est = sg * t$slope_mult_scale * mult_range,
                     lo = sg * t$slope_mult_scale_lo * mult_range, hi = sg * t$slope_mult_scale_hi * mult_range),
          data.frame(method = t$method, metric = met, axis = "add", est = sg * t$slope_add_shift * add_range,
                     lo = sg * t$slope_add_shift_lo * add_range, hi = sg * t$slope_add_shift_hi * add_range))
  }))
  eff$lo2 <- pmin(eff$lo, eff$hi); eff$hi2 <- pmax(eff$lo, eff$hi)
  f <- eff[eff$metric == "Fscore" & eff$axis == "mult", ]; ord <- f$method[order(f$est)]
  ggplot(eff, aes(est, factor(method, levels = ord), colour = method)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lo2, xmax = hi2), width = 0.25, orientation = "y") + geom_point(size = 2) +
    facet_wrap(metric ~ axis, scales = "free_x", ncol = 2,
               labeller = labeller(axis = c(add = sprintf("Additive (0 \u2192 %g)", add_range), mult = sprintf("Multiplicative (0 \u2192 %g)", mult_range)))) +
    labs(x = "Degradation over full severity range", y = NULL) +
    scale_colour_manual(values = c("#00202e", "#003f5c", "#2c4875", "#8a508f", "#bc5090", "#ff6361", "#ff8531", "#ffa600", "#439775")) +
    theme_bw()
}

## Replicate-averaged BEGONE heatmap (one facet per method). Empty tiles mean
## the statistic was undefined in every replicate of that cell.
begone_heatmap <- function(df, value, pal = BEGONE_PALETTE[[value]]) {
  agg <- aggregate(df[[value]], by = list(method = df$method, add_shift = df$add_shift, mult_scale = df$mult_scale), FUN = mean, na.rm = TRUE)
  names(agg)[4] <- "value"
  ggplot(agg, aes(factor(add_shift), factor(mult_scale), fill = value)) +
    geom_tile(linewidth = 0.25) + facet_wrap(~ method) +
    scale_fill_gradient(low = pal[["low"]], high = pal[["high"]], na.value = pal[["na"]]) +
    labs(x = "additive (mean) shift", y = "multiplicative (variance) scale", fill = value, title = paste0("BEGONE \u2014 ", value)) +
    theme_bw(base_size = BASE_SIZE) +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key.width = unit(0.45, "lines"),
          legend.key.height = unit(1.6, "lines"), axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 2.5),
          axis.text.y = element_text(size = BASE_SIZE - 2.5), strip.text = element_text(size = BASE_SIZE - 0.5, margin = margin(2, 2, 2, 2)),
          panel.spacing = unit(3, "pt"), plot.title = element_text(size = BASE_SIZE + 0.5, face = "bold"), plot.margin = margin(4, 5, 2, 3))
}

## Iterations needed before a metric's estimate stabilises: the first run of
## `window` consecutive iterations whose cumulative SD changes by less than
## `tol`, ignoring the first `burn_in` points (small-sample flukes).
convergence_iterations <- function(metric_by_iter, tol = 0.05, window = 8, burn_in = 15) {
  x <- as.numeric(metric_by_iter)
  cum_sd <- sapply(seq_along(x), function(k) sd(x[seq_len(k)]))
  pct_change <- c(NA, abs(diff(cum_sd)) / (abs(head(cum_sd, -1)) + 1e-8))
  below <- !is.na(pct_change) & pct_change < tol
  below[seq_len(min(burn_in, length(below)))] <- FALSE
  first_stable <- NA_integer_
  if (length(below) >= window) {
    run_ok <- sapply(seq_len(length(below) - window + 1), function(i) all(below[i:(i + window - 1)]))
    hit <- which(run_ok)[1]
    if (!is.na(hit)) first_stable <- hit + window - 1   # index of the LAST point in the stable run
  }
  list(cum_sd = cum_sd, pct_change = pct_change, converged_at = first_stable, converged = !is.na(first_stable))
}

## ---------------------------------------------------------------------------
## Paired statistics and BRAID figures
## ---------------------------------------------------------------------------
signif_stars <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns")))))

## Paired Wilcoxon signed-rank test of set_b vs set_a per method (paired on
## iteration), BH-adjusted across the methods in `d`. Non-finite values are
## replaced by `na_as` first (e.g. an F-score with no DE features called).
braid_paired_stats <- function(d, value, set_a = "batch-ve", set_b = "batch+ve", na_as = 0) {
  if (!is.na(na_as)) d[[value]][!is.finite(d[[value]])] <- na_as
  out <- do.call(rbind, lapply(unique(d$method), function(m) {
    s <- d[d$method == m, ]
    a <- s[s$set == set_a, ]; a <- a[order(a$iteration), value]
    b <- s[s$set == set_b, ]; b <- b[order(b$iteration), value]
    n <- min(length(a), length(b)); a <- a[seq_len(n)]; b <- b[seq_len(n)]
    ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]; n <- length(a)
    p <- if (n < 2) NA_real_ else tryCatch(suppressWarnings(wilcox.test(b, a, paired = TRUE))$p.value, error = function(e) NA_real_)
    data.frame(method = m, n_pairs = n, median_diff = if (n) median(b - a) else NA_real_, p_value = p)
  }))
  out$p_adj <- p.adjust(out$p_value, method = "BH"); out
}

## batch-ve vs batch+ve boxplots per method for one metric, split by missingness
## mode, with paired-Wilcoxon brackets. Stars are computed per method and mode
## (bh_across_methods = FALSE reproduces the published figures; TRUE adjusts
## across the methods within each mode instead).
plot_braid_summary <- function(d, metric, title, bh_across_methods = FALSE) {
  lv <- intersect(names(MISS_LABEL), unique(d$missingness))
  d[[metric]][!is.finite(d[[metric]])] <- 0
  d$set <- factor(d$set, levels = c("batch-ve", "batch+ve"))
  d$mvmode <- factor(MISS_LABEL[d$missingness], levels = MISS_LABEL[lv])
  d$combo <- factor(paste0(d$mvmode, ": ", d$set), levels = as.vector(t(outer(MISS_LABEL[lv], levels(d$set), paste, sep = ": "))))
  st <- do.call(rbind, lapply(lv, function(mode) {
    s <- d[d$missingness == mode, ]
    r <- if (bh_across_methods) braid_paired_stats(s, metric)
         else do.call(rbind, lapply(unique(s$method), function(m) braid_paired_stats(s[s$method == m, ], metric)))
    r$mvmode <- MISS_LABEL[[mode]]; r
  }))
  ann <- do.call(rbind, lapply(unique(d$method), function(m) {
    s <- d[d$method == m, ]; hi <- max(s[[metric]], na.rm = TRUE); rng <- diff(range(s[[metric]], na.rm = TRUE))
    if (!is.finite(rng) || rng == 0) rng <- 0.1
    x <- st[st$method == m, ]; lbl <- signif_stars(x$p_adj); ns <- lbl == "ns"
    data.frame(method = m, mvmode = x$mvmode, start = paste0(x$mvmode, ": batch-ve"), end = paste0(x$mvmode, ": batch+ve"),
               y_position = hi + rng * ifelse(ns, 0.07, 0.05), label = lbl, vjust = ifelse(ns, -0.6, 0.2))
  }))
  base <- list(scale_fill_manual(values = BRAID_SET_FILL), scale_colour_manual(values = BRAID_SET_COLOUR), guides(colour = "none"),
               labs(x = NULL, y = METRIC_LABEL[[metric]], fill = "Batch condition", title = title),
               scale_y_continuous(expand = expansion(mult = c(0.04, 0.20))), coord_cartesian(clip = "off"), thin_theme(),
               theme(axis.text.x = element_text(angle = 30, hjust = 1, size = BASE_SIZE - 2)))
  ggplot(d, aes(combo, .data[[metric]], fill = set, colour = set)) + geom_boxplot(outlier.size = 0.4, linewidth = 0.4) +
    facet_wrap(~ method, scales = "free_y") + base +
    suppressWarnings(ggsignif::geom_signif(data = ann, inherit.aes = FALSE, manual = TRUE, tip_length = 0.02,   # ggsignif warns spuriously
                                           aes(xmin = start, xmax = end, annotations = label, y_position = y_position, vjust = vjust)))
}

## Parameter sweep: RMSE per parameter value, batch-ve vs batch+ve, with paired
## Wilcoxon stars BH-adjusted across parameter values.
plot_param_sweep <- function(d, metric = "RMSE") {
  st <- braid_paired_stats(transform(d, method = as.character(param_value)), metric)
  st$param_value <- as.numeric(st$method); st$stars <- signif_stars(st$p_adj)
  d <- d[is.finite(d[[metric]]), , drop = FALSE]
  d$set <- factor(d$set, levels = c("batch-ve", "batch+ve")); lv <- sort(unique(d$param_value))
  d$xf <- factor(d$param_value, levels = lv)
  rng <- diff(range(d[[metric]], na.rm = TRUE)); if (!is.finite(rng) || rng == 0) rng <- 1
  ann <- do.call(rbind, lapply(seq_along(lv), function(i) {
    s <- st$stars[st$param_value == lv[i]]; if (!length(s)) return(NULL)
    data.frame(x = i, xmin = i - 0.1875, xmax = i + 0.1875, y = max(d[[metric]][d$param_value == lv[i]]) + rng * 0.05,
               label = s, vjust = if (s == "ns") -0.55 else 0.25)
  }))
  xl <- PARAM_LABEL[[d$param_name[1]]]
  p <- ggplot(d, aes(xf, .data[[metric]], fill = set, colour = set)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.4, width = 0.7) +
    scale_fill_manual(values = BRAID_SET_FILL) + scale_colour_manual(values = BRAID_SET_COLOUR) + guides(colour = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(x = xl, y = METRIC_LABEL[[metric]], fill = NULL, title = sprintf("%s \u2014 %s vs %s", d$method[1], METRIC_LABEL[[metric]], xl)) +
    thin_theme()
  if (!is.null(ann))
    p <- p + geom_segment(data = ann, inherit.aes = FALSE, aes(x = xmin, xend = xmax, y = y, yend = y), linewidth = 0.3) +
      geom_text(data = ann, inherit.aes = FALSE, aes(x = x, y = y, label = label, vjust = vjust), size = 2.6)
  list(plot = p, stats = st)
}

## ---------------------------------------------------------------------------
## Paired log-ratio forest plots (Figure 3, S7, S8)
## ---------------------------------------------------------------------------
## Per method x dataset x BECA x missingness x metric: r_i = log((batch+ve + c) /
## (batch-ve + c)) paired on iteration; Hodges-Lehmann estimate with exact
## Wilcoxon 95% CI; BH over every row. F-score is sign-flipped so positive = worse.
## c = 1% of the median non-zero logFC error (0 for RMSE and F-score).
FOREST_FILL   <- c("#76C457", "#E45742")
FOREST_COLOUR <- c("#2A7C13", "#972828")
braid_logratio_table <- function(util, min_pairs = 3, ci_level = 0.95) {
  util <- util[util$design == "native", ]
  util$mode <- MISS_LABEL[util$missingness]
  pos <- util$logFC_err[is.finite(util$logFC_err) & util$logFC_err > 0]
  pc <- c(RMSE = 0, logFC_err = if (length(pos)) 0.01 * median(pos) else 0, Fscore = 0)
  cat(sprintf("  logFC error pseudocount: %.6g\n", pc[["logFC_err"]]))
  ci_of <- function(x) {
    n <- length(x); if (n < 2) return(rep(NA_real_, 4))
    if (sd(x) == 0) return(c(x[1], x[1], x[1], 1))
    w <- tryCatch(wilcox.test(x, mu = 0, conf.int = TRUE, conf.level = ci_level, exact = n < 50), error = function(e) NULL,
                  warning = function(w) suppressWarnings(wilcox.test(x, mu = 0, conf.int = TRUE, conf.level = ci_level, exact = FALSE)))
    if (!is.null(w)) return(c(unname(w$estimate), w$conf.int, w$p.value))
    tt <- tryCatch(t.test(x, mu = 0, conf.level = ci_level), error = function(e) NULL)
    if (is.null(tt)) c(mean(x), NA, NA, NA) else c(unname(tt$estimate), tt$conf.int, tt$p.value)
  }
  keys <- unique(util[, c("dataset", "beca", "missingness", "mode", "method")])
  fr <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
    s <- merge(keys[k, ], util)
    do.call(rbind, lapply(names(pc), function(met) {
      mm <- merge(s[s$set == "batch-ve", c("iteration", met)], s[s$set == "batch+ve", c("iteration", met)], by = "iteration", suffixes = c(".ve", ".pve"))
      lo <- mm[[paste0(met, ".ve")]]; hi <- mm[[paste0(met, ".pve")]]
      ok <- is.finite(lo) & is.finite(hi) & (lo + pc[[met]]) > 0 & (hi + pc[[met]]) > 0
      ci <- if (sum(ok) >= min_pairs) ci_of(log((hi[ok] + pc[[met]]) / (lo[ok] + pc[[met]]))) else rep(NA_real_, 4)
      data.frame(keys[k, c("dataset", "beca", "mode", "method")], metric = met, pseudocount = pc[[met]],
                 n_total = nrow(mm), n_pairs = sum(ok), est = ci[1], lo = ci[2], hi = ci[3], p_value = ci[4])
    }))
  }))
  flip <- fr$metric == "Fscore"
  fr$log_ratio_oriented <- ifelse(flip, -fr$est, fr$est)
  fr$lo_oriented <- ifelse(flip, -fr$hi, fr$lo); fr$hi_oriented <- ifelse(flip, -fr$lo, fr$hi)
  fr$p_adj <- p.adjust(fr$p_value, method = "BH"); fr$significant <- !is.na(fr$p_adj) & fr$p_adj < 0.05
  lab <- function(m) if (pc[[m]] > 0) sprintf("%s (+%g)", c(RMSE = "RMSE", logFC_err = "logFC error", Fscore = "F-score")[[m]], pc[[m]])
                     else c(RMSE = "RMSE", logFC_err = "logFC error", Fscore = "F-score")[[m]]
  fr$metric_lab <- factor(vapply(fr$metric, lab, ""), levels = vapply(names(pc), lab, ""))
  fr$mode <- factor(fr$mode); fr$dataset <- factor(fr$dataset); fr$beca <- factor(fr$beca)
  ord <- aggregate(log_ratio_oriented ~ method, fr, mean, na.rm = TRUE)
  fr$method <- factor(fr$method, levels = ord$method[order(ord$log_ratio_oriented)])
  rownames(fr) <- NULL; fr
}
plot_logratio_forest <- function(pf, facets, title, subtitle) {
  ds <- levels(pf$dataset)
  ggplot(pf, aes(log_ratio_oriented, method, colour = dataset, shape = dataset)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo_oriented, xmax = hi_oriented), width = 0, orientation = "y", linewidth = 0.4, position = position_dodge(width = 0.5)) +
    geom_point(aes(fill = interaction(significant, dataset)), size = 1.8, stroke = 0.5, position = position_dodge(width = 0.5)) +
    facet_grid(facets, scales = "free_x") +
    scale_colour_manual(name = "Dataset", values = FOREST_COLOUR[seq_along(ds)], breaks = ds) +
    scale_shape_manual(name = "Dataset", values = c(21, 24)[seq_along(ds)], breaks = ds) +
    scale_fill_manual(values = setNames(c(FOREST_FILL[seq_along(ds)], rep("white", length(ds))),
                                        c(paste0("TRUE.", ds), paste0("FALSE.", ds))), guide = "none") +
    labs(x = "log ratio (batch+ve / batch-ve), positive = worse", y = NULL, fill = NULL, title = title, subtitle = subtitle) +
    theme_bw() + theme(legend.position = "bottom")
}

## ---------------------------------------------------------------------------
## Experiment drivers
## ---------------------------------------------------------------------------
## Manuscript settings and defaults. run_BRAID.R / run_BEGONE.R run one of
## these; the figure scripts run whatever they need through the same drivers,
## so a setting behaves identically either way.
BRAID_ITERATIONS <- 10
BRAID_MV_PROP    <- c(0.3, 0.7)          # total MV fraction, MNAR fraction (3:7 MCAR:MNAR)
BRAID_KSAMP      <- 4                    # KNN-sample K (number of sample neighbours)
BRAID_KFEAT      <- 10                   # KNN-feature K (impute package default)
BRAID_SIM        <- list(n_pseudo = 4, shift = 2.5, scale = 0.25, seed = 42)   # single_batch_sim only
BRAID_SWEEP_SETTING <- c("VP", "combat", "standard", "native")                 # the only setting that sweeps
BRAID_SWEEPS <- list(list(dir = "KNNsample_param_sweep",  method = "KNN-sample",  param = "ksamp", values = c(2, 4, 6, 8, 10)),
                     list(dir = "KNNfeature_param_sweep", method = "KNN-feature", param = "kfeat", values = c(5, 10, 15, 20, 25)),
                     list(dir = "QRILC_param_sweep",      method = "QRILC",       param = "sigma", values = c(0.25, 0.5, 0.75, 1, 1.25)))
BRAID_SETTINGS <- rbind(expand.grid(dataset = c("VP", "HP"), beca = c("combat", "limma"), missingness = c("standard", "batchlinked"),
                                    design = "native", stringsAsFactors = FALSE),
                        expand.grid(dataset = "HP", beca = "combat", missingness = c("standard", "batchlinked"),
                                    design = "single_batch_sim", stringsAsFactors = FALSE))
BEGONE_ITERATIONS  <- 10
BEGONE_MISSINGNESS <- "standard"
BEGONE_ADD_LEVELS  <- seq(0, 4.5, by = 0.5)     # additive (mean) shift
BEGONE_MULT_LEVELS <- seq(0, 0.45, by = 0.05)   # multiplicative (variance) scale
BEGONE_FEATURE_SUBSAMPLE <- 0.35                # fraction of features redrawn per replicate
BEGONE_HP_N_PSEUDO <- 4; BEGONE_HP_PSEUDO_SEED <- 42

## One BRAID setting: results, per-setting boxplots, and (for the reference
## setting only) the three parameter sweeps. Already-finished work is skipped
## unless force = TRUE; an interrupted run always resumes from its checkpoint.
braid_run_setting <- function(dataset = "VP", beca = "combat", missingness = "standard", design = "native",
                              iterations = BRAID_ITERATIONS, methods = MVI_METHODS, force = FALSE) {
  stopifnot(dataset %in% c("VP", "HP"), beca %in% c("combat", "limma"),
            missingness %in% names(MISS_TAG), design %in% c("native", "single_batch_sim"))
  out_dir <- braid_setting_dir(dataset, beca, missingness, design, iterations)
  csv <- file.path(out_dir, "fast_utility.csv"); is_sweep <- identical(c(dataset, beca, missingness, design), BRAID_SWEEP_SETTING)
  done <- file.exists(csv) && (!is_sweep || all(file.exists(vapply(BRAID_SWEEPS, function(s) file.path(out_dir, s$dir, paste0(s$method, "_sweep.csv")), ""))))
  if (done && !force) { cat(sprintf("BRAID | %s | %s | %s | %s  [already done, skipping]\n", dataset, design, beca, missingness)); return(invisible(read.csv(csv))) }
  cat(sprintf("BRAID | %s | %s | BECA = %s | missingness = %s\n  output: %s\n", dataset, design, beca, missingness, out_dir))
  set.seed(42)
  dat <- load_dataset(dataset); ext <- NULL
  if (design == "single_batch_sim") { dat <- do.call(make_single_batch_sim, c(list(dat), BRAID_SIM)); ext <- dat$external_truth }
  methods <- available_methods(methods)
  cat(sprintf("  %d features x %d samples | methods: %s | KNN-sample K = %d (batch-aware reference K = %d)\n",
              nrow(dat$df), ncol(dat$df), paste(methods, collapse = ", "), BRAID_KSAMP, batch_aware_k_sample(dat$batch_f, dat$class_f)))
  if (!file.exists(csv) || force) {
    res <- run_braid(dat$df, dat$batch_f, dat$class_f, BRAID_MV_PROP, iterations, methods, missingness, beca,
                     ksamp = BRAID_KSAMP, kfeat = BRAID_KFEAT, external_truth = ext, checkpoint = file.path(out_dir, "checkpoint.rds"))
    res <- cbind(res, dataset = dataset, beca = beca, missingness = missingness, design = design)
    write.csv(res, csv, row.names = FALSE)
    for (m in names(METRIC_LABEL))
      save_fig(plot_braid_summary(res, m, sprintf("%s \u2014 %s (%s, %s, %s)", METRIC_LABEL[[m]], dataset, beca, MISS_LABEL[[missingness]], design)),
               out_dir, tolower(m), 7.09, 5.0)
  } else res <- read.csv(csv)
  if (is_sweep) for (s in BRAID_SWEEPS) {
    sdir <- file.path(out_dir, s$dir); dir.create(sdir, showWarnings = FALSE); scsv <- file.path(sdir, sprintf("%s_sweep.csv", s$method))
    if (file.exists(scsv) && !force) { cat("  [skip]", s$method, "sweep\n"); next }
    sw <- braid_param_sweep(dat$df, dat$batch_f, dat$class_f, BRAID_MV_PROP, iterations, s$method, s$param, s$values,
                            checkpoint_dir = sdir, missingness = missingness, beca = beca)
    write.csv(sw, scsv, row.names = FALSE)
    save_fig(plot_param_sweep(sw)$plot, sdir, sprintf("%s_param_sweep", gsub("-", "_", s$method)), 3.35, 2.6)
  }
  cat("  done:", out_dir, "\n"); invisible(res)
}
braid_run_all <- function(force = FALSE) invisible(lapply(seq_len(nrow(BRAID_SETTINGS)), function(i)
  do.call(braid_run_setting, c(as.list(BRAID_SETTINGS[i, ]), list(force = force)))))

## One BEGONE severity sweep: grid results plus the three heatmaps.
begone_run_dataset <- function(dataset = "VP", iterations = BEGONE_ITERATIONS, methods = MVI_METHODS, force = FALSE) {
  stopifnot(dataset %in% c("VP", "HP"))
  out_dir <- begone_run_dir(dataset, iterations, "combat", BEGONE_MISSINGNESS, length(methods))
  csv <- file.path(out_dir, "begone_grid.csv")
  if (file.exists(csv) && !force) { cat(sprintf("BEGONE | %s  [already done, skipping]\n", dataset)); return(invisible(read.csv(csv))) }
  cat(sprintf("BEGONE | %s | grid %d x %d | output: %s\n", dataset, length(BEGONE_ADD_LEVELS), length(BEGONE_MULT_LEVELS), out_dir))
  set.seed(42)
  dat <- load_dataset(dataset)
  G_ref <- create_ground_truth(dat$df, dat$batch_f, beca = "combat")     # batch-free reference
  batch_f <- dat$batch_f
  if (dataset == "HP") {                                                # relabel AFTER correction
    stopifnot(BEGONE_HP_N_PSEUDO <= min(table(dat$class_f)))
    batch_f <- assign_pseudo_batches(dat$class_f, BEGONE_HP_N_PSEUDO, BEGONE_HP_PSEUDO_SEED)
  }
  print(table(class = dat$class_f, batch = batch_f))
  ksamp <- BRAID_KSAMP# batch_aware_k_sample(batch_f, dat$class_f)
  methods <- available_methods(methods)
  cat(sprintf("  %d features x %d samples | KNN-sample K = %d | methods: %s\n", nrow(G_ref), ncol(G_ref), ksamp, paste(methods, collapse = ", ")))
  grid <- run_begone(G_ref, batch_f, dat$class_f, BEGONE_ADD_LEVELS, BEGONE_MULT_LEVELS, BRAID_MV_PROP, iterations, methods,
                     BEGONE_MISSINGNESS, beca = "combat", ksamp = ksamp, feature_subsample = BEGONE_FEATURE_SUBSAMPLE,
                     checkpoint = file.path(out_dir, "begone_checkpoint.rds"))
  write.csv(grid, csv, row.names = FALSE)
  for (v in names(BEGONE_PALETTE)) save_fig(begone_heatmap(grid, v), out_dir, paste0("begone_heatmap_", v), 7.09, 5.0)
  cat("  done:", out_dir, "\n"); invisible(grid)
}
begone_run_all <- function(force = FALSE) invisible(lapply(c("VP", "HP"), begone_run_dataset, force = force))
