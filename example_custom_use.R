###############################################################################
## example_custom_use.R — a worked example of R/BRAID_custom.R
##
## Runs end to end on a small synthetic dataset, so it needs no input files and
## finishes in a minute or two. Walks through the five things a new user does:
##
##   1. load BRAID and point it at a scratch results folder
##   2. inspect a dataset before running it
##   3. run BRAID on your own matrix
##   4. add your own imputation method and compare it against the built-ins
##   5. add your own batch-effect correction method, and run a BEGONE sweep
##
## Usage, from the repository root:
##   Rscript examples/example_custom_use.R
## In RStudio, just source the file.
##
## Results go to results_example/ and are deleted at the end unless you set
## KEEP_RESULTS <- TRUE below. Nothing here touches results/.
###############################################################################
KEEP_RESULTS <- FALSE

## --------------------------------------------------------------- 1. setup --
## Source the library, then the helpers. braid_setup() does both in one call
## once BRAID_custom.R is loaded, which is handier interactively.
.here <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/", mustWork = TRUE))))
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(dirname(normalizePath(f$ofile, winslash = "/", mustWork = TRUE))))
  getwd()
}
root <- .here()
Sys.setenv(BRAID_RESULTS = file.path(root, "results_example"))
BRAID_ROOT <- root
source(file.path(root, "R", "BRAID_functions.R"))     # the library first
source(file.path(root, "R", "BRAID_custom.R"))        # then the helpers
## Equivalently, in an interactive session, this one call does all three:
##   source("R/BRAID_custom.R"); braid_setup("/path/to/BRAID", results = "...")
cat("BRAID loaded from", root, "\n  results ->", results_path(), "\n")

## ------------------------------------------------- 2. a synthetic dataset --
## Features x samples, log2 scale, complete. 24 samples: 4 batches x 2 classes
## x 3 replicates. 60 of the 300 features are differentially expressed, and
## affine batch effects are injected so there is something for BRAID to correct.
set.seed(2026)
n_feat <- 300; n_samp <- 24
batch <- factor(rep(1:4, each = 6))
class <- factor(rep(rep(c("ctrl", "treat"), each = 3), times = 4))
base  <- rnorm(n_feat, mean = 20, sd = 2)
mat   <- matrix(rnorm(n_feat * n_samp, mean = base, sd = 0.5), n_feat, n_samp,
                dimnames = list(sprintf("feature_%03d", seq_len(n_feat)), sprintf("sample_%02d", seq_len(n_samp))))
de <- 1:60                                                     # the features that truly differ
mat[de, class == "treat"] <- mat[de, class == "treat"] + 1.5
mat <- sim_batch_effects(mat, batch, shift = 1.5, scale = 0.15, seed = 1)

## ------------------------------------------------------ 3. look before you run
## braid_check() reports completeness, scale, group sizes, the K it would
## choose, and which F-score truth applies. Feature names carry no spike-in
## tags here, so the truth is the DE set called on the corrected ground truth.
cat("\n== braid_check ==\n")
chk <- braid_check(mat, batch, class)

## ----------------------------------------------------------- 4. run BRAID --
## Three built-in methods, three iterations: enough to see the shape of the
## output. Drop the `methods` argument to run all nine.
cat("\n== BRAID, built-in methods ==\n")
res <- braid_on_data(mat, batch, class,
                     methods = c("KNN-sample", "Mean", "SVD"),
                     iterations = 3)
cat("\nMean over iterations:\n")
print(aggregate(cbind(RMSE, logFC_err, Fscore) ~ method + set, res, mean, na.action = na.pass), row.names = FALSE)

## ------------------------------------------- 5. add your own imputer --------
## The function receives a features x samples matrix containing NAs and returns
## it filled. Keep the `...`: BRAID passes ksamp, kfeat, ntrees and sigma to
## every method, and a function without `...` will error when they arrive.
## This one fills each missing value with its feature's median.
median_impute <- function(mat, ...) {
  mat <- as.matrix(mat)
  med <- apply(mat, 1, median, na.rm = TRUE); med[is.na(med)] <- mean(mat, na.rm = TRUE)
  i <- which(is.na(mat)); mat[i] <- med[row(mat)[i]]
  mat
}
register_imputer("Median", median_impute)
print(list_imputers()$custom)

cat("\n== BRAID, with the custom method alongside a built-in ==\n")
res2 <- braid_on_data(mat, batch, class,
                      methods = c("Mean", "Median"),
                      iterations = 3,
                      out_dir = results_path("example_braid"))     # writes CSV + figures
cat("\nMean over iterations:\n")
print(aggregate(cbind(RMSE, logFC_err, Fscore) ~ method + set, res2, mean, na.action = na.pass), row.names = FALSE)

## -------------------------------------- 6. add your own correction method --
## fn(mat, batch, class) returns the corrected matrix. The registered name then
## works anywhere BRAID takes `beca`. This one centres each batch on the
## overall feature mean — a crude stand-in for ComBat.
register_corrector("meancentre", function(mat, batch, class = NULL) {
  mat <- as.matrix(mat); gm <- rowMeans(mat, na.rm = TRUE)
  for (b in unique(batch)) {
    i <- batch == b
    mat[, i] <- mat[, i] - rowMeans(mat[, i, drop = FALSE], na.rm = TRUE) + gm
  }
  mat
})
cat("\n== BRAID, using the custom corrector as the BECA ==\n")
res3 <- braid_on_data(mat, batch, class, methods = "Mean", iterations = 2, beca = "meancentre")
print(aggregate(RMSE ~ method + set, res3, mean), row.names = FALSE)

## ---------------------------------------------------------- 7. BEGONE -------
## BEGONE injects known batch effects, so it needs a batch-free matrix: correct
## the data first. The grid here is deliberately tiny (3 x 3 cells, 2
## replicates); the manuscript uses 10 x 10 cells and 10 replicates.
cat("\n== BEGONE sweep on a small grid ==\n")
clean <- apply_beca(mat, batch, "combat", class)
grid <- begone_on_data(clean, batch, class,
                       methods = c("Mean", "Median"),
                       iterations = 2,
                       add_levels  = c(0, 2, 4),
                       mult_levels = c(0, 0.2, 0.4),
                       out_dir = results_path("example_begone"))
cat("\nMean batch+ve RMSE by severity:\n")
print(aggregate(batch_pve ~ method + add_shift + mult_scale, grid, mean), row.names = FALSE)

## ---------------------------------------------------------------- cleanup --
cat("\nDone.\n")
if (!KEEP_RESULTS) {
  unlink(file.path(root, "results_example"), recursive = TRUE)
  cat("Removed results_example/. Set KEEP_RESULTS <- TRUE at the top to keep it.\n")
} else cat("Output kept in", file.path(root, "results_example"), "\n")
