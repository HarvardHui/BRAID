###############################################################################
## S05_begone.R — Supplementary Figures S21-S28 and Tables S1-S2 (BEGONE)
##
## Runs the BEGONE sweeps for VP and HP itself (skipping either if it is
## already finished), then draws the figures. Nothing has to be edited.
##   Figures S21-S23: VP heatmaps of replicate-averaged batch+ve RMSE, logFC error, F-score
##   Figures S24-S26: the same for HP
##   Figure S27:      HP per-method degradation over the full severity range
##   Figure S28:      HP variance decomposition
##   Table S1:        VP runs with a missing (non-finite) metric, per method
##   Table S2:        the same per severity cell
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
begone_run_all()
grids <- lapply(c(VP = "VP", HP = "HP"), function(ds) read.csv(file.path(begone_run_dir(ds), "begone_grid.csv")))

## ------------------------------------------------ Figures S21-S26: heatmaps
ids <- list(VP = c(batch_pve = "S21", logfc_pve = "S22", fscore_pve = "S23"), HP = c(batch_pve = "S24", logfc_pve = "S25", fscore_pve = "S26"))
for (ds in names(grids)) for (v in names(BEGONE_PALETTE))
  save_fig(begone_heatmap(grids[[ds]], v), fig_dir, sprintf("Fig%s_begone_heatmap_%s_%s", ids[[ds]][[v]], v, ds), 7.09, 5.0)

## ------------------------------------------- Figures S27, S28: HP statistics
q <- begone_quant(grids$HP)
write.csv(do.call(rbind, q$variance), file.path(tab_dir, "begone_HP_variance_decomposition.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(q$slopes), function(o) cbind(metric = BEGONE_OUTCOMES[[o]], q$slopes[[o]]$table))),
          file.path(tab_dir, "begone_HP_slopes.csv"), row.names = FALSE)
save_fig(plot_begone_degradation(q$slopes, max(grids$HP$add_shift), max(grids$HP$mult_scale)), fig_dir, "FigS27_begone_degradation_HP", 9, 6)
save_fig(plot_variance_decomposition(q$variance), fig_dir, "FigS28_begone_variance_decomposition_HP", 6, 4)

## ------------------------------------ Tables S1, S2: incomplete VP metrics
g <- grids$VP; metrics <- names(BEGONE_OUTCOMES)
n_expected <- length(unique(g$method)) * length(unique(g$add_shift)) * length(unique(g$mult_scale)) * length(unique(g$iteration))
cat(sprintf("VP grid: %d of %d expected rows present\n", nrow(g), n_expected))
t1 <- aggregate(g[metrics], by = list(method = g$method), FUN = function(x) sum(!is.finite(x)))
t1 <- t1[rowSums(t1[metrics]) > 0, ]; names(t1)[-1] <- BEGONE_OUTCOMES
t2 <- do.call(rbind, lapply(metrics, function(m) {
  bad <- g[!is.finite(g[[m]]), ]; if (!nrow(bad)) return(NULL)
  a <- aggregate(list(n_runs = rep(1, nrow(bad))), by = bad[c("method", "add_shift", "mult_scale")], FUN = sum)
  cbind(metric = BEGONE_OUTCOMES[[m]], a[order(-a$n_runs), ])
}))
write.csv(t1, file.path(tab_dir, "TableS1_begone_VP_missing_metrics_by_method.csv"), row.names = FALSE)
write.csv(t2, file.path(tab_dir, "TableS2_begone_VP_missing_metrics_by_cell.csv"), row.names = FALSE)
print(t1, row.names = FALSE)
cat("Figures S21-S28 written to", fig_dir, "; Tables S1-S2 to", tab_dir, "\n")
