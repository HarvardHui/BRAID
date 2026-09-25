###############################################################################
## S03_mice_rf_limitations.R — Supplementary Figures S4, S5 and S6
##
## MICE-PMM, MICE-norm and RF are run with FEATURES as variables (conventional;
## mice sees t(mat), with quickpred(mincor = 0.4, minpuc = 0.3) to keep it
## tractable) and with SAMPLES as variables (used in the manuscript), on a
## random 200-feature subset of VP and HP, with batch effects present (raw data)
## and removed (ComBat with the class design). Missing values: 30% total, 3:7
## MCAR:MNAR, replicate i seeded with 123456 + i; 10 replicates per cell.
## mice: m = 1, maxit = 5; RF: missForest, ntree = 100.
##   Figure S4: fraction of missing cells filled.
##   Figure S5: wall-clock seconds per imputation.
##   Figure S6: RMSE over the cells actually filled (point size = fraction filled).
## Output: results/figures/supplementary/, results/mice_rf_limitations/
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
fig_dir <- results_path("figures", "supplementary"); out_dir <- results_path("mice_rf_limitations")

DATASETS <- c(VP = "Van Puyvelde", HP = "Haslett & Pescatori")
N_REPS <- 10; SEED0 <- 123456; MV_PROP <- c(0.3, 0.7); N_FEATURES <- 200; FEATURE_SEED <- 42
METHODS <- c("MICE-PMM", "MICE-norm", "RF"); MICE_ENGINE <- c("MICE-PMM" = "pmm", "MICE-norm" = "norm")
VARIANTS <- c("features-as-variables" = TRUE, "samples-as-variables" = FALSE)   # value = transpose (and use quickpred)
FILL   <- c(No = BRAID_SET_FILL[["batch-ve"]],   Yes = BRAID_SET_FILL[["batch+ve"]])
COLOUR <- c(No = BRAID_SET_COLOUR[["batch-ve"]], Yes = BRAID_SET_COLOUR[["batch+ve"]])
thin <- theme_bw(base_size = 9) + theme(legend.position = "right", legend.key.size = unit(0.8, "lines"), panel.grid.minor = element_blank(),
                                        axis.text.x = element_text(angle = 15, hjust = 1), plot.subtitle = element_text(size = 8, colour = "grey30"))

## One imputation: time, cells filled, and RMSE (all MV cells if complete; filled cells otherwise).
run_one <- function(mat_true, transpose, method, seed, batch_f) {
  mis <- simulate_mvs(mat_true, MV_PROP, batch_factor = batch_f, seed = seed)$msdata
  work <- if (transpose) t(mis) else mis
  keep <- colSums(!is.na(work)) > 0
  sub <- as.data.frame(work[, keep, drop = FALSE])
  colnames(sub) <- make.unique(make.names(tolower(gsub("[^A-Za-z0-9]+", "_", colnames(sub)))))
  t0 <- Sys.time()
  out <- tryCatch(if (method == "RF") as.matrix(missForest::missForest(sub, ntree = 100)$ximp) else {
    args <- list(data = sub, m = 1, maxit = 5, method = MICE_ENGINE[[method]], printFlag = FALSE)
    if (transpose) args$pred <- tryCatch(suppressWarnings(mice::quickpred(sub, mincor = 0.4, minpuc = 0.3)), error = function(e) NULL)
    as.matrix(mice::complete(do.call(mice::mice, args)))
  }, error = function(e) NULL)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  imp <- NULL
  if (!is.null(out)) { filled <- work; filled[, keep] <- out; imp <- if (transpose) t(filled) else filled; dimnames(imp) <- dimnames(mis) }
  was <- is.na(mis); got <- was & if (is.null(imp)) FALSE else !is.na(imp)
  rmse_on <- function(mask) if (is.null(imp) || !any(mask)) NA_real_ else { m <- mis; m[!mask] <- 0; m[mask] <- NA; Rmse(imp, m, mat_true) }
  data.frame(seconds = secs, n_missing = sum(was), n_filled = sum(got), frac_filled = sum(got) / sum(was),
             rmse = if (sum(got) == sum(was)) rmse_on(was) else NA_real_, rmse_filled = rmse_on(got))
}

rows <- list()
for (ds in names(DATASETS)) {
  d <- load_dataset(ds); d$df <- as.matrix(d$df)
  if (nrow(d$df) > N_FEATURES) { set.seed(FEATURE_SEED); d$df <- d$df[sort(sample(nrow(d$df), N_FEATURES)), , drop = FALSE] }
  corrected <- tryCatch(as.matrix(sva::ComBat(d$df, batch = d$batch_f, mod = model.matrix(~ d$class_f))),
                        error = function(e) as.matrix(sva::ComBat(d$df, batch = d$batch_f)))
  for (mth in METHODS) for (vn in names(VARIANTS)) for (bs in c("Yes", "No")) {
    for (i in seq_len(N_REPS))
      rows[[length(rows) + 1L]] <- cbind(data.frame(dataset = DATASETS[[ds]], method = mth, variant = vn, batch = bs, replicate = i),
                                         run_one(if (bs == "Yes") d$df else corrected, VARIANTS[[vn]], mth, SEED0 + i, d$batch_f))
    cat(sprintf("%s | %s | %s | batch = %s done\n", ds, mth, vn, bs))
  }
}
res <- do.call(rbind, rows)
res$batch <- factor(res$batch, levels = c("No", "Yes")); res$method <- factor(res$method, levels = METHODS)
res$variant <- factor(res$variant, levels = names(VARIANTS))
write.csv(res, file.path(out_dir, "mice_rf_limitations_raw.csv"), row.names = FALSE)
write.csv(aggregate(cbind(seconds, frac_filled, rmse, rmse_filled) ~ dataset + method + variant + batch, res, mean, na.action = na.pass),
          file.path(out_dir, "mice_rf_time_and_rmse.csv"), row.names = FALSE)

## Figure S4 — fraction of missing cells filled
p4 <- ggplot(res, aes(dataset, frac_filled, fill = batch, colour = batch)) +
  geom_boxplot(aes(group = interaction(dataset, batch)), width = 0.5, outlier.shape = NA, position = position_dodge(width = 0.6), linewidth = 0.4, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6), alpha = 0.8, size = 1.3) +
  scale_fill_manual(values = FILL) + scale_colour_manual(values = COLOUR) + guides(colour = "none") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), limits = c(0, NA)) + facet_grid(method ~ variant) +
  labs(x = NULL, y = "Missing cells filled", fill = "Batch effect\npresent?",
       title = "How much of the missingness was resolved \u2014 MICE-PMM, MICE-norm and RF",
       subtitle = "100% = fully imputed; anything below is a partial imputation") + thin
save_fig(p4, fig_dir, "FigS4_mice_rf_fraction_filled", 7.09, 6.4)

## Figure S5 — run time
p5 <- ggplot(res, aes(dataset, seconds, fill = batch, colour = batch)) + geom_boxplot(outlier.size = 0.5, linewidth = 0.4) +
  scale_fill_manual(values = FILL) + scale_colour_manual(values = COLOUR) + guides(colour = "none") +
  facet_grid(method ~ variant, scales = "free_y") + scale_y_log10() +
  labs(x = NULL, y = "Seconds per imputation", fill = "Batch effect\npresent?",
       title = "Imputation run time \u2014 MICE-PMM, MICE-norm and RF", subtitle = "one point per replicate; log-scaled axis") + thin
save_fig(p5, fig_dir, "FigS5_mice_rf_time", 7.09, 6.4)

## Figure S6 — RMSE over the cells actually filled
okp <- res[res$n_filled > 0 & is.finite(res$rmse_filled), ]
p6 <- ggplot(okp, aes(dataset, rmse_filled, colour = batch)) +
  geom_point(aes(size = frac_filled), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.7), alpha = 0.75) +
  scale_colour_manual(values = COLOUR) + scale_size_continuous(range = c(0.5, 3), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1)) +
  facet_grid(method ~ variant, scales = "free_y") +
  labs(x = NULL, y = "RMSE (cells actually filled)", colour = "Batch effect\npresent?", size = "Fraction\nfilled",
       title = "Imputation error, including partial imputations \u2014 MICE-PMM, MICE-norm and RF",
       subtitle = "small points filled few cells, so their error is not comparable to a complete run") + thin
save_fig(p6, fig_dir, "FigS6_mice_rf_rmse_filled", 7.09, 6.4)
cat("Figures S4-S6 written to", fig_dir, "\n")
