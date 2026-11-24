###############################################################################
## BRAID_custom.R — using BRAID on your own data, and with your own methods
##
## An optional add-on. Source it AFTER R/BRAID_functions.R; it adds helpers and
## a method registry without changing anything the manuscript scripts rely on,
## so results from run_BRAID.R, run_BEGONE.R and the figure scripts are
## unaffected whether or not this file is loaded.
##
##   source("R/BRAID_functions.R")       # BRAID_ROOT must be set first; see below
##   source("R/BRAID_custom.R")
##
## WHAT IT ADDS
##   braid_setup()        point BRAID at a repository and load everything
##   braid_check()        inspect a matrix before running: completeness,
##                        negatives, group sizes, which DE truth will be used
##   braid_on_data()      run BRAID on your own matrix, with sane defaults
##   begone_on_data()     run a BEGONE severity sweep on your own matrix
##   register_imputer()   add an imputation method without editing the library
##   register_corrector() add a batch-effect correction method likewise
##   list_imputers() / list_correctors() / unregister_imputer()
##
## MINIMAL EXAMPLE
##   source("R/BRAID_functions.R"); source("R/BRAID_custom.R")
##   braid_check(mat, batch, class)                    # look before you leap
##   res <- braid_on_data(mat, batch, class, iterations = 10)
##
## ADDING YOUR OWN IMPUTER
##   register_imputer("MyMethod", function(mat, ...) my_impute(mat))
##   res <- braid_on_data(mat, batch, class, methods = c("KNN-sample", "MyMethod"))
##
## The function receives a features x samples matrix containing NAs and must
## return a matrix of the same shape with the NAs filled. Keep the `...`: BRAID
## passes its hyperparameters (ksamp, kfeat, ntrees, sigma) to every method.
##
## REQUIREMENTS FOR YOUR DATA
##   * features x samples, log-transformed, and COMPLETE (no NAs). BRAID
##     defines its ground truth by correcting the complete matrix and then
##     simulating missingness into it, so any pre-existing NAs must be removed
##     first, e.g. mat <- mat[complete.cases(mat), ].
##   * batch and class are vectors of length ncol(mat).
##   * Every batch should contain every class, or batch and biology are
##     confounded and the correction cannot separate them. braid_check() warns.
###############################################################################

## ---------------------------------------------------------------------------
## Setup: locate the repository and load the library
## ---------------------------------------------------------------------------
## braid_setup("/path/to/BRAID") sets BRAID_ROOT, sources the library and this
## file, and optionally redirects where results are written.
braid_setup <- function(root, results = NULL) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "R", "BRAID_functions.R")))
    stop("No R/BRAID_functions.R under ", root)
  assign("BRAID_ROOT", root, envir = globalenv())
  if (!is.null(results)) Sys.setenv(BRAID_RESULTS = results)
  source(file.path(root, "R", "BRAID_functions.R"), local = globalenv())
  source(file.path(root, "R", "BRAID_custom.R"), local = globalenv())
  .braid_install_wrappers()
  cat("BRAID loaded from", root, "\n  results ->", results_path(), "\n")
  invisible(root)
}

## ---------------------------------------------------------------------------
## Method registries
## ---------------------------------------------------------------------------
if (!exists(".braid_registry", inherits = FALSE))
  .braid_registry <- list(imputers = list(), correctors = list())

## Add an imputer. fn(mat, ...) takes a features x samples matrix with NAs and
## returns the filled matrix. Names registered here take precedence over the
## built-ins, so an existing name needs overwrite = TRUE.
register_imputer <- function(name, fn, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, is.function(fn))
  if (name %in% MVI_METHODS && !overwrite)
    stop("'", name, "' is a built-in method; pass overwrite = TRUE to shadow it.")
  if (!is.null(.braid_registry$imputers[[name]]) && !overwrite)
    stop("'", name, "' is already registered; pass overwrite = TRUE to replace it.")
  .braid_registry$imputers[[name]] <<- fn
  cat("registered imputer:", name, "\n"); invisible(name)
}
unregister_imputer <- function(name) { .braid_registry$imputers[[name]] <<- NULL; invisible(name) }
list_imputers <- function() list(builtin = MVI_METHODS, custom = names(.braid_registry$imputers))

## Add a batch-effect correction method. fn(mat, batch, class) returns the
## corrected matrix. Use the registered name wherever BRAID takes `beca`.
register_corrector <- function(name, fn, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, is.function(fn))
  if (name %in% c("combat", "limma") && !overwrite)
    stop("'", name, "' is built in; pass overwrite = TRUE to shadow it.")
  .braid_registry$correctors[[name]] <<- fn
  cat("registered corrector:", name, "\n"); invisible(name)
}
list_correctors <- function() list(builtin = c("combat", "limma"), custom = names(.braid_registry$correctors))

## Wrap the library's dispatchers so registered methods are tried first.
## Installed by .braid_install_wrappers(), which runs at the bottom of this file
## and again from braid_setup(). It is safe to call repeatedly: the originals
## are captured once, so wrappers never nest, and it does nothing (returning
## FALSE) when R/BRAID_functions.R has not been loaded yet.
.braid_install_wrappers <- function(env = globalenv()) {
  if (!exists("impute_one", envir = env) || !exists("apply_beca", envir = env)) return(invisible(FALSE))
  if (!exists(".braid_impute_one_orig", envir = env, inherits = FALSE))
    assign(".braid_impute_one_orig", get("impute_one", envir = env), envir = env)
  if (!exists(".braid_apply_beca_orig", envir = env, inherits = FALSE))
    assign(".braid_apply_beca_orig", get("apply_beca", envir = env), envir = env)

  assign("impute_one", function(mat, method, ...) {
    fn <- .braid_registry$imputers[[method]]
    if (is.null(fn)) return(.braid_impute_one_orig(mat, method, ...))
    mat <- as.matrix(mat)
    res <- tryCatch(fn(mat, ...),
                    error = function(e) { warning("Could not apply ", method, "; leaving it unimputed. Reason: ",
                                                  conditionMessage(e)); mat })
    res <- as.matrix(res)
    if (!identical(dim(res), dim(mat)))
      stop("Imputer '", method, "' returned ", nrow(res), "x", ncol(res), "; expected ", nrow(mat), "x", ncol(mat), ".")
    dimnames(res) <- dimnames(mat); res
  }, envir = env)

  assign("apply_beca", function(data, batch, method = c("combat", "limma"), class = NULL) {
    m <- if (is.character(method)) method[1] else method
    fn <- if (is.character(m)) .braid_registry$correctors[[m]] else NULL
    if (is.null(fn)) return(.braid_apply_beca_orig(data, batch, m, class))
    out <- as.matrix(fn(as.matrix(data), batch, class))
    dimnames(out) <- dimnames(as.matrix(data)); out
  }, envir = env)
  invisible(TRUE)
}

## Methods known to BRAID right now: built-ins whose packages are installed,
## plus everything registered here.
available_methods_custom <- function(methods) {
  custom <- intersect(methods, names(.braid_registry$imputers))
  c(available_methods(setdiff(methods, custom)), custom)
}

## ---------------------------------------------------------------------------
## Inspect a dataset before running
## ---------------------------------------------------------------------------
## Reports what BRAID will do with the matrix, and warns about the conditions
## that most often cause confusing results. Returns its findings invisibly.
braid_check <- function(mat, batch, class, quiet = FALSE) {
  mat <- as.matrix(mat); batch <- as.factor(batch); class <- as.factor(class)
  say <- function(...) if (!quiet) cat(...)
  ok <- TRUE
  say(sprintf("Matrix: %d features x %d samples\n", nrow(mat), ncol(mat)))
  if (length(batch) != ncol(mat) || length(class) != ncol(mat)) {
    say("  ERROR: batch and class must have one entry per column.\n"); ok <- FALSE
  }
  n_na <- sum(is.na(mat))
  if (n_na) { say(sprintf("  ERROR: %d missing values (%.2f%%). BRAID needs a complete matrix; try mat[complete.cases(mat), ].\n",
                          n_na, 100 * n_na / length(mat))); ok <- FALSE }
  if (is.null(rownames(mat))) { say("  ERROR: rownames are required (scoring aligns methods on feature names).\n"); ok <- FALSE }
  rng <- range(mat, na.rm = TRUE)
  say(sprintf("  Range: %.2f to %.2f%s\n", rng[1], rng[2],
              if (rng[2] > 1000) "  <- looks unlogged; BRAID expects log-transformed data" else ""))
  if (rng[1] < 0) say(sprintf("  Note: %d negative values. Fine for most methods, but NMF-based ones need non-negative input.\n",
                              sum(mat < 0)))
  tab <- table(batch, class)
  say("  Samples per batch x class:\n"); if (!quiet) print(tab)
  if (any(tab == 0)) { say("  WARNING: some batch x class combinations are empty; batch and biology are confounded.\n"); ok <- FALSE }
  k <- max(1L, min(tab) - 1L)
  say(sprintf("  Smallest batch-class group: %d  (batch-aware KNN-sample K would be %d)\n", min(tab), k))
  ts <- .de_truth_source(rownames(mat))
  say(sprintf("  F-score truth: %s\n", if (ts == "spikein")
    "spike-in (YEAST/ECOLI differential, HUMAN not), from feature names" else
    "reference (DE called on the batch-corrected ground truth)"))
  say(if (ok) "  Ready.\n" else "  Fix the items above before running.\n")
  invisible(list(ok = ok, n_missing = n_na, range = rng, table = tab, ksamp_batch_aware = k, de_truth = ts))
}

## ---------------------------------------------------------------------------
## Run BRAID / BEGONE on your own data
## ---------------------------------------------------------------------------
## Same engine as the manuscript runs, with arguments that make sense for a new
## dataset. Returns the results table; writes nothing unless out_dir is given.
braid_on_data <- function(mat, batch, class, methods = MVI_METHODS, iterations = 10,
                          mv_prop = c(0.3, 0.7), missingness = "standard", beca = "combat",
                          ksamp = NULL, kfeat = 10, seed0 = 123456, out_dir = NULL, check = TRUE, ...) {
  mat <- as.matrix(mat); batch <- as.factor(batch); class <- as.factor(class)
  chk <- braid_check(mat, batch, class, quiet = !check)
  if (check && !chk$ok) stop("braid_check() found problems; fix them or pass check = FALSE to proceed anyway.")
  if (is.null(ksamp)) { ksamp <- chk$ksamp_batch_aware; cat("Using batch-aware KNN-sample K =", ksamp, "\n") }
  methods <- available_methods_custom(methods)
  cat("Methods:", paste(methods, collapse = ", "), "\n")
  res <- run_braid(mat, batch, class, mv_prop, iterations, methods, missingness, beca,
                   ksamp = ksamp, kfeat = kfeat, seed0 = seed0,
                   checkpoint = if (is.null(out_dir)) NULL else file.path(out_dir, "checkpoint.rds"), ...)
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(res, file.path(out_dir, "fast_utility.csv"), row.names = FALSE)
    for (m in names(METRIC_LABEL))
      save_fig(plot_braid_summary(cbind(res, missingness = missingness), m, METRIC_LABEL[[m]]), out_dir, tolower(m), 7.09, 5.0)
    cat("Written to", out_dir, "\n")
  }
  res
}

## clean_mat must be batch-free (BEGONE injects known batch effects into it).
## Correct your data first if it has native batch structure:
##   clean <- apply_beca(mat, batch, "combat", class)
begone_on_data <- function(clean_mat, batch, class, methods = MVI_METHODS, iterations = 10,
                           add_levels = seq(0, 4.5, by = 0.5), mult_levels = seq(0, 0.45, by = 0.05),
                           mv_prop = c(0.3, 0.7), missingness = "standard", beca = "combat",
                           ksamp = NULL, feature_subsample = 0.35, seed0 = 424242, out_dir = NULL, check = TRUE, ...) {
  clean_mat <- as.matrix(clean_mat); batch <- as.factor(batch); class <- as.factor(class)
  chk <- braid_check(clean_mat, batch, class, quiet = !check)
  if (check && !chk$ok) stop("braid_check() found problems; fix them or pass check = FALSE to proceed anyway.")
  if (is.null(ksamp)) ksamp <- chk$ksamp_batch_aware
  methods <- available_methods_custom(methods)
  cat(sprintf("BEGONE grid: %d x %d cells x %d replicates | KNN-sample K = %d\n",
              length(add_levels), length(mult_levels), iterations, ksamp))
  grid <- run_begone(clean_mat, batch, class, add_levels, mult_levels, mv_prop, iterations, methods,
                     missingness, beca, ksamp = ksamp, feature_subsample = feature_subsample, seed0 = seed0,
                     checkpoint = if (is.null(out_dir)) NULL else file.path(out_dir, "begone_checkpoint.rds"), ...)
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(grid, file.path(out_dir, "begone_grid.csv"), row.names = FALSE)
    for (v in names(BEGONE_PALETTE)) save_fig(begone_heatmap(grid, v), out_dir, paste0("begone_heatmap_", v), 7.09, 5.0)
    cat("Written to", out_dir, "\n")
  }
  grid
}

if (.braid_install_wrappers())
  cat("BRAID_custom.R loaded. See braid_setup(), braid_check(), braid_on_data(), register_imputer().\n") else
  cat("BRAID_custom.R loaded before R/BRAID_functions.R; call braid_setup(\"<repo>\") to finish loading.\n")
