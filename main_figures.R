###############################################################################
## main_figures.R — every main-manuscript figure, in one script.
##
## Run this and nothing else: it runs the BRAID settings and BEGONE sweeps the
## figures need (skipping any that are already finished), then draws Figures
## 2-6. Nothing in this file has to be edited. Use run_BRAID.R / run_BEGONE.R
## only if you want one setting on its own.
##
## The limma BRAID settings are run because the Figure 3 q-values are
## BH-adjusted jointly over every setting. Figures 2 and 6B are computed here;
## Figure 2 needs no input data at all.
##
## Output: results/figures/main/ (figures), results/tables/ (numbers behind them)
###############################################################################
## Locate the repository: the script's own folder (Rscript, source() or RStudio),
## walking up to the folder that holds R/BRAID_functions.R. Override with BRAID_ROOT.
# .script_dir <- function() {
#   a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
#   if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = TRUE)))
#   for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile, winslash = "/", mustWork = TRUE)))
#   if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
#     p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
#     if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/", mustWork = TRUE)))
#   }
#   getwd()
# }

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
fig_dir <- results_path("figures", "main"); tab_dir <- results_path("tables")

## ---------------------------------------------------------------------------
## Prerequisites: 8 native BRAID settings (Figure 3) and BEGONE VP (Figures 4, 5)
## ---------------------------------------------------------------------------
native <- BRAID_SETTINGS[BRAID_SETTINGS$design == "native", ]
invisible(lapply(seq_len(nrow(native)), function(i) do.call(braid_run_setting, as.list(native[i, ]))))
begone_run_dataset("VP")

## ---------------------------------------------------------------------------
## Figure 2 - what an injected batch effect looks like. A Gaussian matrix
## (1000 features x 20 samples, 5 batches of 4) receives affine batch effects at
## the four combinations of additive shift 0 / 2.5 and multiplicative scale
## 0 / 0.25 (0, 0 is the unaffected control). Each combination gives a PCA of the
## samples and the per-batch distribution of expression values.
## ---------------------------------------------------------------------------
FIG2_ROWS <- 1000; FIG2_COLS <- 20; FIG2_BATCHES <- 5
FIG2_DATA_SEED <- 1234567890; FIG2_SIM_SEED <- 1
FIG2_ADD <- c(0, 2.5); FIG2_MULT <- c(0, 0.25)
FIG2_COLS_PAL <- c("#FF9999", "#99CCFF", "#99FFB3", "#FFE699", "#CC99FF")
set.seed(FIG2_DATA_SEED)
row_means <- rnorm(FIG2_ROWS, mean = 10, sd = 1); row_sds <- pmax(rnorm(FIG2_ROWS, mean = 1, sd = 0.2), 0.1)
fig2_mat <- matrix(nrow = FIG2_ROWS, ncol = FIG2_COLS)
for (i in 1:FIG2_ROWS) fig2_mat[i, ] <- rnorm(FIG2_COLS, mean = row_means[i], sd = row_sds[i])
dimnames(fig2_mat) <- list(paste0("Gene_", 1:FIG2_ROWS), paste0("Sample_", 1:FIG2_COLS))
fig2_bf <- rep(seq_len(FIG2_BATCHES), each = FIG2_COLS / FIG2_BATCHES)
fig2_theme <- theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(size = 10))
for (loc in FIG2_ADD) for (mul in FIG2_MULT) {
  m <- sim_batch_effects(fig2_mat, fig2_bf, shift = loc, scale = mul, seed = FIG2_SIM_SEED)
  ttl <- sprintf("Additive = %s, multiplicative = %s", loc, mul)
  tag <- sprintf("add%s_mult%s", loc, mul)
  pca <- summary(prcomp(t(m), scale. = TRUE, center = TRUE))
  pc <- as.data.frame(pca$x); pc$batch <- as.factor(fig2_bf)
  save_fig(ggplot(pc, aes(PC1, PC2, fill = batch)) +
             geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
             geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
             geom_point(size = 3, shape = 21, color = "black") + scale_fill_manual(values = FIG2_COLS_PAL) +
             labs(x = sprintf("PC1 (%s%%)", pca$importance[2, 1] * 100), y = sprintf("PC2 (%s%%)", pca$importance[2, 2] * 100), title = ttl) +
             fig2_theme, fig_dir, paste0("Fig2_pca_", tag), 4.2, 3.0)
  den <- data.frame(value = as.numeric(m), batch = as.factor(rep(fig2_bf, each = nrow(m))))
  save_fig(ggplot(den, aes(value, fill = batch)) + geom_density(alpha = 0.3) +
             scale_fill_manual(values = FIG2_COLS_PAL) + scale_y_continuous(expand = expansion(mult = 0)) +
             labs(x = "Expression value", y = "Density", title = ttl) + fig2_theme,
           fig_dir, paste0("Fig2_density_", tag), 4.2, 3.0)
}

## ---------------------------------------------------------------------------
## Figure 3 — paired log ratio (batch+ve / batch-ve) of RMSE, logFC error and
## F-score; VP and HP; sample-wise and BEAMs missingness; ComBat.
## ---------------------------------------------------------------------------
fr <- braid_logratio_table(load_braid_results())
write.csv(fr, file.path(tab_dir, "logratio_summary.csv"), row.names = FALSE)
pf <- fr[is.finite(fr$log_ratio_oriented), ]
forest_sub <- paste0("filled points are significant after BH correction (95% CI, wilcoxon)",
                     if (any(fr$pseudocount > 0)) "\n(+c) marks a pseudocount: ratio of (metric + c), shrunk toward zero" else "")
save_fig(plot_logratio_forest(pf[pf$beca == "combat", ], mode ~ metric_lab, "Relative degradation under batch effects, by setting", forest_sub),
         fig_dir, "Fig3_logratio_forest_combat", 7.09, 7.0)

## ---------------------------------------------------------------------------
## Figure 4 — variance decomposition of the VP BEGONE grid.
## Figure 5A — per-method degradation over the full severity range (VP).
## Figure 5B — replicate-averaged logFC error heatmap, all nine methods (VP).
## ---------------------------------------------------------------------------
grid <- read.csv(file.path(begone_run_dir("VP"), "begone_grid.csv"))
q <- begone_quant(grid)
write.csv(do.call(rbind, q$variance), file.path(tab_dir, "begone_VP_variance_decomposition.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(q$slopes), function(o) cbind(metric = BEGONE_OUTCOMES[[o]], q$slopes[[o]]$table))),
          file.path(tab_dir, "begone_VP_slopes.csv"), row.names = FALSE)
save_fig(plot_variance_decomposition(q$variance), fig_dir, "Fig4_begone_variance_decomposition_VP", 6, 4)
save_fig(plot_begone_degradation(q$slopes, max(grid$add_shift), max(grid$mult_scale)), fig_dir, "Fig5A_begone_degradation_VP", 9, 6)
save_fig(begone_heatmap(grid, "logfc_pve"), fig_dir, "Fig5B_begone_heatmap_logfc_VP", 7.09, 5.0)

## ---------------------------------------------------------------------------
## Figure 6A — KNN-sample RMSE against K, batch-ve vs batch+ve (VP, ComBat,
## sample-wise missingness), paired Wilcoxon, BH across K.
## ---------------------------------------------------------------------------
sw <- plot_param_sweep(read.csv(file.path(braid_setting_dir("VP", "combat", "standard"), "KNNsample_param_sweep", "KNN-sample_sweep.csv")))
write.csv(sw$stats, file.path(tab_dir, "param_sweep_KNN-sample_rmse_signif.csv"), row.names = FALSE)
save_fig(sw$plot, fig_dir, "Fig6A_param_sweep_KNN-sample", 3.35, 2.6)

## ---------------------------------------------------------------------------
## Figure 6B — % of the K neighbours used by KNN-sample that come from another
## batch, at increasing K and MV proportion (VP). Neighbours are chosen exactly
## as do.knn() chooses them (euclidean distance between samples, candidates =
## samples observed for that feature); the ordering is computed once per
## missing value and read off at every K. batch-ve (ComBat-corrected) is the
## negative control; dashed lines give chance level for the candidate pools.
## ---------------------------------------------------------------------------
K_VALUES <- c(2, 4, 6, 8, 10, 15, 20); MV_PROPS <- c(0.1, 0.2, 0.3, 0.4, 0.5); MNAR_FRAC <- 0.7
XB_ITER <- 10; XB_SEED0 <- 42; FEATURE_CAP <- 2000
cross_batch_rates <- function(mis, bnum, ks) {
  kmax <- max(ks); dist_mat <- as.matrix(dist(t(mis), method = "euclidean")); na_mat <- is.na(mis)
  cross_sum <- numeric(length(ks)); n_cells <- 0L; chance_sum <- 0
  for (i in which(rowSums(na_mat) > 0)) {
    obs <- which(!na_mat[i, ]); mvs <- which(na_mat[i, ])
    if (length(obs) < kmax) next                       # do.knn() errors when the pool is smaller than K
    ord <- matrix(t(apply(dist_mat[mvs, obs, drop = FALSE], 1, order)), nrow = length(mvs))
    cross <- matrix(bnum[obs][ord], nrow = length(mvs)) != bnum[mvs]
    cum <- matrix(t(apply(cross[, seq_len(kmax), drop = FALSE], 1, cumsum)), nrow = length(mvs))
    cross_sum <- cross_sum + colSums(cum[, ks, drop = FALSE])
    chance_sum <- chance_sum + sum(vapply(mvs, function(s) mean(bnum[obs] != bnum[s]), numeric(1)))
    n_cells <- n_cells + length(mvs)
  }
  if (!n_cells) return(NULL)
  list(rate = 100 * (cross_sum / ks) / n_cells, chance = 100 * chance_sum / n_cells)
}
vp <- load_VP(); bnum <- as.integer(vp$batch_f)
mats <- list("batch+ve" = as.matrix(vp$df), "batch-ve" = as.matrix(create_ground_truth(vp$df, vp$batch_f, beca = "combat")))
xb <- list()
for (path in names(mats)) for (p in MV_PROPS) for (it in seq_len(XB_ITER)) {
  seed <- XB_SEED0 + it
  mis <- filter_high_mv_features(simulate_mvs(mats[[path]], c(p, MNAR_FRAC), seed = seed)$msdata, 0.6)
  if (nrow(mis) > FEATURE_CAP) { set.seed(seed); mis <- mis[sort(sample(nrow(mis), FEATURE_CAP)), , drop = FALSE] }
  r <- cross_batch_rates(mis, bnum, K_VALUES)
  if (!is.null(r)) xb[[length(xb) + 1L]] <- data.frame(path = path, mv_prop = p, iteration = it, K = K_VALUES, pct_cross = r$rate, chance = r$chance)
}
xb <- do.call(rbind, xb)
agg <- merge(aggregate(cbind(pct_cross, chance) ~ path + mv_prop + K, xb, mean),
             setNames(aggregate(pct_cross ~ path + mv_prop + K, xb, sd), c("path", "mv_prop", "K", "sd_cross")))
write.csv(agg, file.path(tab_dir, "cross_batch_neighbours_summary.csv"), row.names = FALSE)
agg$mv_lab <- factor(sprintf("%d%%", round(100 * agg$mv_prop)), levels = sprintf("%d%%", round(100 * MV_PROPS)))
agg$path <- factor(agg$path, levels = c("batch+ve", "batch-ve"))
p6b <- ggplot(agg, aes(K, pct_cross, colour = mv_lab, group = mv_lab)) +
  geom_hline(data = aggregate(chance ~ path, agg, mean), aes(yintercept = chance), linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pct_cross - sd_cross, ymax = pct_cross + sd_cross), width = 0.35, linewidth = 0.3, alpha = 0.65) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.4) + facet_wrap(~ path) +
  scale_x_continuous(breaks = K_VALUES) + scale_colour_viridis_d(end = 0.88, option = "D") +
  labs(x = "K (number of sample neighbours)", y = "% of K neighbours from another batch", colour = "Missing values",
       title = "KNN-sample draws its neighbours from within batch",
       subtitle = paste("dashed line = chance level given the candidate pools; error bars = SD over", XB_ITER, "MV masks")) +
  theme_bw(base_size = 9) + theme(legend.position = "right", legend.key.size = unit(0.8, "lines"), panel.grid.minor = element_blank(),
                                  plot.subtitle = element_text(size = 8, colour = "grey30"))
save_fig(p6b, fig_dir, "Fig6B_knn_cross_batch_neighbours", 6.5, 3.8)
cat("Main figures written to", fig_dir, "\n")
