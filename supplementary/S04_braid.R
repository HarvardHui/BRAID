###############################################################################
## S04_braid.R — Supplementary Figures S7-S20 (BRAID)
##
## Runs every BRAID setting of the manuscript itself (skipping any already
## finished), then draws the figures. Nothing has to be edited.
##   native:           VP, HP x combat, limma x standard, batchlinked
##   single_batch_sim: HP x combat x standard, batchlinked
## plus the parameter sweeps, which belong to VP / combat / standard / native.
##   Figure S7:     paired log-ratio forest plot, limma
##   Figure S8:     F-score log ratios without NMFBatch, ComBat and limma
##   Figures S9-S19: batch-ve vs batch+ve boxplots with paired Wilcoxon tests
##   Figure S20:    KNN-feature K (A) and QRILC sigma (B) sweeps
## Output: results/figures/supplementary/, results/tables/
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
fig_dir <- results_path("figures", "supplementary"); tab_dir <- results_path("tables")
braid_run_all()
util <- load_braid_results()

## ------------------------------------------------- Figures S7, S8: forests
fr <- braid_logratio_table(util)
pf <- fr[is.finite(fr$log_ratio_oriented), ]
save_fig(plot_logratio_forest(pf[pf$beca == "limma", ], mode ~ metric_lab, "Relative degradation under batch effects, by setting",
                              paste0("filled points are significant after BH correction (95% CI, wilcoxon)",
                                     if (any(fr$pseudocount > 0)) "\n(+c) marks a pseudocount: ratio of (metric + c), shrunk toward zero" else "")),
         fig_dir, "FigS7_logratio_forest_limma", 7.09, 7.0)
pfs <- pf[pf$metric == "Fscore" & pf$method != "NMFBatch", ]; pfs$method <- droplevels(pfs$method)
save_fig(plot_logratio_forest(pfs, mode ~ beca, "F-score: relative degradation under batch effects",
                              "excludes NMFBatch; filled points are significant after BH correction (95% CI, wilcoxon)"),
         fig_dir, "FigS8_logratio_forest_fscore_no_nmfbatch", 7.09, 7.0)

## ------------------------------------------- Figures S9-S19: boxplots
## Each figure pools the two missingness modes of one dataset x BECA x design.
## Settings without a supplementary figure number are still drawn ("extra").
fig_id <- c("HP_combat_single_batch_sim_RMSE" = "S9",  "HP_combat_single_batch_sim_logFC_err" = "S10",
            "VP_combat_native_RMSE" = "S11", "HP_combat_native_RMSE" = "S12", "VP_combat_native_logFC_err" = "S13",
            "HP_combat_native_logFC_err" = "S14", "VP_combat_native_Fscore" = "S15", "HP_combat_native_Fscore" = "S16",
            "VP_limma_native_RMSE" = "S17", "VP_limma_native_logFC_err" = "S18", "VP_limma_native_Fscore" = "S19")
groups <- unique(util[, c("dataset", "beca", "design")])
for (g in seq_len(nrow(groups))) {
  d <- merge(groups[g, ], util)
  for (m in names(METRIC_LABEL)) {
    key <- paste(groups$dataset[g], groups$beca[g], groups$design[g], m, sep = "_")
    name <- sprintf("Fig%s_braid_%s", if (key %in% names(fig_id)) fig_id[[key]] else "extra", key)
    title <- sprintf("%s \u2014 %s (%s, %s)", METRIC_LABEL[[m]], groups$dataset[g], groups$beca[g], if (groups$design[g] == "native") "nosim" else "sim")
    save_fig(plot_braid_summary(d, m, title), fig_dir, name, 7.09, 5.0)
    st <- do.call(rbind, lapply(unique(d$missingness), function(mode) cbind(missingness = mode, braid_paired_stats(d[d$missingness == mode, ], m))))
    write.csv(st, file.path(tab_dir, paste0("braid_signif_", key, ".csv")), row.names = FALSE)
  }
}

## ------------------------------------------------- Figure S20: sweeps
for (k in 1:2) {
  s <- list(c("KNNfeature_param_sweep", "KNN-feature"), c("QRILC_param_sweep", "QRILC"))[[k]]
  sw <- plot_param_sweep(read.csv(file.path(braid_setting_dir("VP", "combat", "standard"), s[1], paste0(s[2], "_sweep.csv"))))
  write.csv(sw$stats, file.path(tab_dir, sprintf("param_sweep_%s_rmse_signif.csv", s[2])), row.names = FALSE)
  save_fig(sw$plot, fig_dir, sprintf("FigS20%s_param_sweep_%s", LETTERS[k], s[2]), 3.35, 2.6)
}
cat("Figures S7-S20 written to", fig_dir, "\n")
