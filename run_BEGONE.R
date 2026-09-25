###############################################################################
## run_BEGONE.R - run ONE BEGONE severity sweep (manuscript Sections 2.2, 3.2)
##
## For readers who want a single dataset. To reproduce every main-manuscript
## figure in one go, run main_figures.R instead; it runs whatever it needs.
##
## Usage, from anywhere:
##   Rscript run_BEGONE.R [DATASET]        DATASET = VP | HP (default VP)
## In RStudio, just source the file (edit the default below to change dataset).
##
## The native batch structure is first removed with ComBat to give a batch-free
## reference. Affine batch effects are then injected over a grid of additive
## shifts (0-4.5) x multiplicative scales (0-0.45), 10 replicates per cell, each
## on a fresh 35% subsample of features. VP keeps its 4 instrument batches; HP
## (2 studies) is relabelled into 4 class-stratified pseudo-batches after
## correction. KNN-sample K is batch-aware: smallest batch-class group - 1.
##
## Output: results/begone/begone_10iter_<DATASET>_combat_standard_9methods/
##   begone_grid.csv    one row per method x cell x replicate (RMSE, logFC
##                      error and F-score for batch-ve and batch+ve)
##   begone_heatmap_{batch_pve,logfc_pve,fscore_pve} (.png/.pdf)
## Settings live in R/BRAID_functions.R (BEGONE_ADD_LEVELS, BEGONE_MULT_LEVELS, ...).
## An interrupted sweep resumes from its begone_checkpoint.rds.
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
DATASET <- if (length(args) >= 1) args[1] else "VP"

begone_run_dataset(DATASET, force = TRUE)
