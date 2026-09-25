###############################################################################
## S01_batch_effect_strength.R — Supplementary Figures S1 and S2
##
## A Gaussian matrix (2000 features x 40 samples, 4 batches of 10) receives
## affine batch effects along one axis at a time: additive shift 0-4.5 or
## multiplicative scale 0-0.45. Replicate i uses seed 1000 + i at every level
## (common random numbers), 100 replicates per level. Batch-effect strength is
## quantified on the top 10 PCs by batch silhouette, variance-weighted PC
## regression R^2, and the variance fraction on PC1.
##   Figure S2: mean strength over the first 10 replicates, 95% t-interval.
##   Figure S1: replicates needed before the paired difference between each
##              pair of adjacent levels is, and stays, significant (95% CI
##              of the mean difference excludes zero).
## Needs no input data. Output: results/figures/supplementary/, results/batch_effect_strength/
###############################################################################
## Locate the repository: the script's own folder (Rscript, source() or RStudio),
## walking up to the folder that holds R/BRAID_functions.R. Override with BRAID_ROOT.
.script_dir <- function() {
  # 1. Rscript: --file=/path/to/script.R
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/",
      mustWork = TRUE
    )))
  }
  
  # 2. RStudio / source(): look for the sourced file
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    f <- frames[[i]]$ofile
    if (!is.null(f)) {
      return(dirname(normalizePath(
        f,
        winslash = "/",
        mustWork = TRUE
      )))
    }
  }
  
  # 3. RStudio: currently active source document
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    ctx <- tryCatch(
      rstudioapi::getActiveDocumentContext(),
      error = function(e) NULL
    )
    
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(
        ctx$path,
        winslash = "/",
        mustWork = TRUE
      )))
    }
  }
  
  # 4. Last resort
  getwd()
}

BRAID_ROOT <- Sys.getenv("BRAID_ROOT", unset = "")
if (!nzchar(BRAID_ROOT)) { d <- .script_dir(); for (i in 1:3) { if (file.exists(file.path(d, "R", "BRAID_functions.R"))) break; d <- dirname(d) }; BRAID_ROOT <- d }
if (!file.exists(file.path(BRAID_ROOT, "R", "BRAID_functions.R")))
  stop("Cannot find R/BRAID_functions.R. Run the script from the repository, or set BRAID_ROOT.")
source(file.path(BRAID_ROOT, "R", "BRAID_functions.R"))
fig_dir <- results_path("figures", "supplementary"); out_dir <- results_path("batch_effect_strength")

SEED_DATA <- 1234567890; N_ROWS <- 2000; N_COLS <- 40; N_BATCH <- 4
SWEEPS <- list(additive = list(param = "shift", levels = seq(0, 4.5, by = 0.5)),
               multiplicative = list(param = "scale", levels = seq(0, 0.45, by = 0.05)))
N_ITER <- 10; N_ITER_MAX <- 100; SEED_BASE <- 1000; N_PCS <- 10; CI_LEVEL <- 0.95
METRIC_NAMES <- c(batch_sil = "Batch silhouette", pcr_weighted = "PC regression (variance-weighted R\u00b2)", pc1_var_frac = "Variance on PC1")
METRICS <- names(METRIC_NAMES)
strength_theme <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 8, colour = "grey30"))

## ---------------------------------------------------------- strength measures
batch_silhouette <- function(pcs, batch) {
  batch <- as.factor(batch); d <- as.matrix(dist(pcs))
  mean(vapply(seq_len(nrow(pcs)), function(i) {
    own <- which(batch == batch[i] & seq_len(nrow(pcs)) != i); if (!length(own)) return(0)
    a <- mean(d[i, own])
    b <- min(vapply(setdiff(levels(batch), as.character(batch[i])), function(g) mean(d[i, which(batch == g)]), numeric(1)))
    (b - a) / max(a, b)
  }, numeric(1)))
}
strength_metrics <- function(mat, batch) {
  X <- t(mat); X <- X[, apply(X, 2, function(v) is.finite(sd(v)) && sd(v) > 0), drop = FALSE]
  pca <- prcomp(X, center = TRUE, scale. = TRUE); k <- min(N_PCS, ncol(pca$x))
  pcs <- pca$x[, seq_len(k), drop = FALSE]; vp <- (pca$sdev^2)[seq_len(k)]
  r2 <- vapply(seq_len(k), function(j) if (sd(pcs[, j]) == 0) 0 else summary(lm(pcs[, j] ~ as.factor(batch)))$r.squared, numeric(1))
  data.frame(batch_sil = batch_silhouette(pcs, batch), pcr_weighted = sum(vp / sum(vp) * r2), pc1_var_frac = vp[1] / sum(pca$sdev^2))
}

## ------------------------------------------------------ baseline and sweeps
set.seed(SEED_DATA)
row_means <- rnorm(N_ROWS, mean = 10, sd = 1); row_sds <- pmax(rnorm(N_ROWS, mean = 1, sd = 0.2), 0.1)
newmat <- matrix(rnorm(N_ROWS * N_COLS, mean = row_means, sd = row_sds), nrow = N_ROWS, ncol = N_COLS,
                 dimnames = list(paste0("Gene_", seq_len(N_ROWS)), paste0("Sample_", seq_len(N_COLS))))
bf <- rep(seq_len(N_BATCH), each = N_COLS / N_BATCH)
raw <- list()
for (sw in names(SWEEPS)) {
  p <- SWEEPS[[sw]]$param
  for (lv in SWEEPS[[sw]]$levels) for (i in seq_len(N_ITER_MAX)) {
    sim <- sim_batch_effects(newmat, bf, shift = if (p == "shift") lv else 0, scale = if (p == "scale") lv else 0, seed = SEED_BASE + i)
    raw[[length(raw) + 1L]] <- cbind(data.frame(param = p, level = lv, iteration = i), strength_metrics(sim, bf))
  }
  cat(sw, "sweep complete\n")
}
raw <- do.call(rbind, raw)
write.csv(raw, file.path(out_dir, "strength_raw_iterations.csv"), row.names = FALSE)

## ------------------------------------------ Figure S2: strength vs severity
tcrit <- function(n) qt(1 - (1 - CI_LEVEL) / 2, df = n - 1)
r10 <- raw[raw$iteration <= N_ITER, ]
summ <- do.call(rbind, lapply(split(r10, list(r10$param, r10$level), drop = TRUE), function(s) do.call(rbind, lapply(METRICS, function(m) {
  x <- s[[m]][is.finite(s[[m]])]; n <- length(x); hw <- tcrit(n) * sd(x) / sqrt(n)
  data.frame(param = s$param[1], level = s$level[1], metric = m, n = n, mean = mean(x), sd = sd(x), ci_lower = mean(x) - hw, ci_upper = mean(x) + hw)
}))))
write.csv(summ, file.path(out_dir, "strength_summary_ci.csv"), row.names = FALSE)
summ$metric_lab <- factor(METRIC_NAMES[summ$metric], levels = unname(METRIC_NAMES))
for (k in seq_along(SWEEPS)) {
  sw <- names(SWEEPS)[k]
  p <- ggplot(summ[summ$param == SWEEPS[[sw]]$param, ], aes(level, mean)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#894b77", alpha = 0.25) +
    geom_line(colour = "#894b77", linewidth = 0.8) + geom_point(colour = "#894b77", size = 1.6) +
    facet_wrap(~ metric_lab, scales = "free_y") +
    labs(x = sprintf("%s batch-effect strength (%s)", sw, SWEEPS[[sw]]$param), y = "Quantified batch-effect strength",
         title = sprintf("Batch-effect strength vs injected %s severity", sw),
         subtitle = sprintf("mean of %d paired replicates; band = %.0f%% CI", N_ITER, 100 * CI_LEVEL)) + strength_theme
  save_fig(p, fig_dir, sprintf("FigS2%s_strength_vs_%s", LETTERS[k], sw), 7.09, 3)
}

## ------------------------- Figure S1: replicates to separate adjacent levels
first_sustained <- function(ok) {   # first index from which `ok` stays TRUE to the end
  ok[is.na(ok)] <- FALSE; if (!any(ok)) return(NA_integer_)
  which(rev(cumprod(rev(as.integer(ok)))) == 1L)[1]
}
excludes_zero <- function(x) if (sd(x) == 0) mean(x) != 0 else abs(mean(x)) > tcrit(length(x)) * sd(x) / sqrt(length(x))
sep <- list()
for (p in unique(raw$param)) {
  lv <- sort(unique(raw$level[raw$param == p]))
  for (m in METRICS) for (k in seq_len(length(lv) - 1L)) {
    mm <- merge(raw[raw$param == p & raw$level == lv[k], c("iteration", m)], raw[raw$param == p & raw$level == lv[k + 1], c("iteration", m)],
                by = "iteration", suffixes = c(".lo", ".hi"))
    d <- mm[[paste0(m, ".hi")]] - mm[[paste0(m, ".lo")]]; d <- d[is.finite(d)]
    if (length(d) < 3) next
    i <- first_sustained(vapply(2:length(d), function(n) isTRUE(excludes_zero(d[seq_len(n)])), logical(1)))
    dn <- d[seq_len(min(N_ITER, length(d)))]
    sep[[length(sep) + 1L]] <- data.frame(param = p, metric = m, level_lo = lv[k], level_hi = lv[k + 1], level_mid = (lv[k] + lv[k + 1]) / 2,
                                          mean_diff = mean(dn), separable_at_n10 = isTRUE(excludes_zero(dn)),
                                          n_required = if (is.na(i)) NA_integer_ else (2:length(d))[i])
  }
}
sep <- do.call(rbind, sep)
write.csv(sep, file.path(out_dir, "adjacent_level_separability.csv"), row.names = FALSE)
sep$metric_lab <- factor(METRIC_NAMES[sep$metric], levels = unname(METRIC_NAMES))
for (k in seq_along(SWEEPS)) {
  sw <- names(SWEEPS)[k]; d <- sep[sep$param == SWEEPS[[sw]]$param & is.finite(sep$n_required), ]
  p <- ggplot(d, aes(level_mid, n_required)) +
    geom_hline(yintercept = N_ITER, linetype = "dotted", colour = "grey30", linewidth = 0.4) +
    geom_line(colour = "#894b77", linewidth = 0.5) + geom_point(colour = "#894b77", size = 1.6) +
    facet_wrap(~ metric_lab, scales = "free_y") +
    labs(x = sprintf("%s severity (midpoint of the adjacent pair)", sw), y = "Replicates needed to separate the pair",
         title = sprintf("Cost of resolving adjacent %s severity levels", sw),
         subtitle = sprintf("dotted line = n = %d; pairs above it are unresolved at that count", N_ITER)) + strength_theme
  save_fig(p, fig_dir, sprintf("FigS1%s_separability_cost_%s", LETTERS[k], sw), 7.09, 3)
}
cat("Figures S1 and S2 written to", fig_dir, "\n")
