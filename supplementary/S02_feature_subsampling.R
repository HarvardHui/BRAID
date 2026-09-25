###############################################################################
## S02_feature_subsampling.R — Supplementary Figure S3
##
## How small can the per-iteration feature subsample be before the RMSE
## estimate becomes noise, and how many iterations does it then need? On the
## Van Puyvelde dataset, KNN-sample (K = 4) is run on fractions of 0.05, 0.10,
## 0.20, 0.35, 0.50 and 1.00 of the features, 30 replicates each, against the
## ComBat-corrected ground truth with 30% missing values (3:7 MCAR:MNAR).
## Replicate r uses seed 90000 + r. Because VP carries spike-in tags, the
## subsample is stratified to keep YEAST/ECOLI proteins represented; datasets
## without such tags fall back to plain random sampling.
## Convergence = the first run of 8 consecutive iterations whose cumulative SD
## changes by less than 5%, ignoring the first 15.
##   Figure S3A: cumulative SD of RMSE against iteration, per fraction
##   Figure S3B: running mean of RMSE against iteration, per fraction
## Needs only the VP input data.
## Output: results/figures/supplementary/, results/feature_subsampling/
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

fig_dir <- results_path("figures", "supplementary"); out_dir <- results_path("feature_subsampling")

FRACTIONS <- c(0.05, 0.10, 0.20, 0.35, 0.50, 1.00)
N_REPS    <- 30                 # iterations swept per fraction
METHOD    <- "KNN-sample"       # fast, low-variance method for this diagnostic
MV_PROP   <- c(0.3, 0.7); KSAMP <- 4
TOL       <- 0.05               # "converged" = <5% relative change in cumulative SD
SEED0     <- 1234567890; REP_SEED0 <- 90000

## Stratified on the spike-in pattern where present, plain random otherwise.
sample_features_safe <- function(mat, n, spike_pattern = "YEAST|ECOLI") {
  if (n >= nrow(mat)) return(mat)
  is_diff <- grepl(spike_pattern, rownames(mat), ignore.case = TRUE)
  idx <- if (any(is_diff)) {
    keep_diff <- sample(which(is_diff), min(sum(is_diff), round(n * 0.35)))
    sort(c(keep_diff, sample(which(!is_diff), min(sum(!is_diff), n - length(keep_diff)))))
  } else sort(sample(nrow(mat), n))
  mat[idx, , drop = FALSE]
}

set.seed(SEED0)
dat <- load_dataset("VP")
raw <- as.matrix(dat$df); batch_factor <- as.factor(dat$batch_f); class_factor <- as.factor(dat$class_f)
cat(if (any(grepl("YEAST|ECOLI", rownames(raw), ignore.case = TRUE)))
      "Spike-in tags (YEAST/ECOLI) detected \u2014 using species-stratified subsampling.\n"
    else "No spike-in tags detected \u2014 using plain random feature subsampling.\n")
GT <- create_ground_truth(raw, batch_factor, beca = "combat", class = class_factor)

run_one <- function(frac, rep_seed) {
  G_sub <- sample_features_safe(GT, max(50, round(nrow(GT) * frac)))
  mv    <- simulate_mvs_mode(G_sub, MV_PROP, batch_factor, "standard", seed = rep_seed)
  misA  <- filter_high_mv_features(mv$msdata, 0.6)
  imp   <- impute_methods(misA, METHOD, verbose = FALSE, ksamp = KSAMP)[[METHOD]]
  Rmse(imp, misA, G_sub[rownames(misA), , drop = FALSE])
}

results <- NULL
for (frac in FRACTIONS) {
  n_feat <- max(50, round(nrow(GT) * frac))
  cat(sprintf("\nFraction = %.2f (~%d features)\n", frac, n_feat))
  vals <- rep(NA_real_, N_REPS)
  for (r in seq_len(N_REPS)) {
    vals[r] <- tryCatch(run_one(frac, REP_SEED0 + r), error = function(e) NA_real_)
    if (r %% 10 == 0) cat("  rep", r, "of", N_REPS, "\n")
  }
  vals <- vals[!is.na(vals)]
  if (length(vals) < 2) { message("  too many failed reps at this fraction \u2014 skipping"); next }
  conv <- convergence_iterations(vals, tol = TOL)
  results <- rbind(results, data.frame(fraction = frac, n_features = n_feat, iteration = seq_along(vals),
                                       running_mean = cumsum(vals) / seq_along(vals), cum_sd = conv$cum_sd,
                                       pct_change = conv$pct_change, converged_at = conv$converged_at))
}
summary_tbl <- unique(results[, c("fraction", "n_features", "converged_at")])
cat("\n== Iterations needed to converge (tol = 5%), by feature fraction ==\n"); print(summary_tbl, row.names = FALSE)
write.csv(results, file.path(out_dir, "feature_subsample_robustness.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(out_dir, "feature_subsample_summary.csv"), row.names = FALSE)

sub_plot <- function(y, ylab, title)   # cum_sd is undefined at iteration 1, so drop non-finite points
  ggplot(results[is.finite(results[[y]]), ], aes(iteration, .data[[y]], colour = factor(fraction))) + geom_line() +
    labs(x = "Iteration", y = ylab, colour = "Feature\nfraction", title = title) + theme_minimal(base_size = 12)
save_fig(sub_plot("cum_sd", "Cumulative SD of RMSE", "Cumulative SD vs. iteration count, by feature fraction"),
         fig_dir, "FigS3A_feature_subsample_cum_sd", 8, 5)
save_fig(sub_plot("running_mean", "Running mean RMSE", sprintf("%s RMSE: running mean vs. iteration count, by feature fraction", METHOD)),
         fig_dir, "FigS3B_feature_subsample_running_mean", 8, 5)
cat("Figure S3 written to", fig_dir, "\n")
