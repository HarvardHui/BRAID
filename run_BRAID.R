###############################################################################
## run_BRAID.R - run ONE BRAID setting (manuscript Sections 2.1, 3.1, 3.3)
##
## For readers who want a single setting. To reproduce every main-manuscript
## figure in one go, run main_figures.R instead; it runs whatever it needs.
##
## Usage, from anywhere:
##   Rscript run_BRAID.R [DATASET] [BECA] [MISSINGNESS] [DESIGN]
##     DATASET      VP | HP                        (default VP)
##     BECA         combat | limma                 (default combat)
##     MISSINGNESS  standard | batchlinked         (default standard; batchlinked = BEAMs)
##     DESIGN       native | single_batch_sim      (default native)
## In RStudio, just source the file (edit the defaults below to change setting).
##
## native:            ground truth = BECA applied to the complete multi-batch data.
## single_batch_sim:  the largest native batch is split into 4 class-stratified
##                    pseudo-batches, affine batch effects (shift 2.5, scale 0.25)
##                    are injected, and the pre-injection data are the ground truth.
##
## Output: results/braid/braid_[sim_]10iter_<DATASET>_<BECA>_<standardmissing|beams>/
##   fast_utility.csv   one row per method x batch condition x iteration
##   rmse, logfc_err, fscore (.png/.pdf)   batch-ve vs batch+ve boxplots
##   *_param_sweep/     KNN-sample K, KNN-feature K and QRILC sigma sweeps; run
##                      only for VP / combat / standard / native (Section 3.3)
## Settings live in R/BRAID_functions.R (BRAID_ITERATIONS, BRAID_KSAMP, ...).
## An interrupted run resumes from its checkpoint.rds; a finished one is redone.
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

args <- commandArgs(trailingOnly = TRUE)
DATASET     <- if (length(args) >= 1) args[1] else "VP"
BECA        <- if (length(args) >= 2) args[2] else "combat"
MISSINGNESS <- if (length(args) >= 3) args[3] else "standard"
DESIGN      <- if (length(args) >= 4) args[4] else "native"

braid_run_setting(DATASET, BECA, MISSINGNESS, DESIGN, force = TRUE)
