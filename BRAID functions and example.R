# Load required libraries
library(proBatch)
library(ggfortify)
library(ggplot2)
library(corrplot)
library(gtools)
library(mice)
library(sva)
library(vegan)
library(imputeLCMD)
library(limma)
library(reshape2)
library(missForest)
library(ggrepel)
library(dplyr)
library(grid)
library(gridExtra)
library(pROC)
library(ggpubr)
library(magrittr)
library(ggplotify)
library(pheatmap)
library(tidyr)
library(here)


wd <- here()


source(paste0(wd,"/BRAID_BEGONE/backend functions.R"))

# ========== BRAID: Modular Workflow for MVI Evaluation ==========

# --------------------- Step 1: Initialization ---------------------
# These functions are specificially for the Van Puyvelde dataset used in this example
# Load complete expression dataset with batch effects
load_design_data <- function(wd, files) {
  lapply(files, function(f) read.table(paste0(wd, f), sep = "\t", header = TRUE))
}

load_intensity_data <- function(wd, files) {
  lapply(files, function(f) read.table(paste0(wd, f), sep = "\t", header = TRUE))
}

extract_exprs <- function(intensity_data) {
  intensity_data[, c(1, 3:ncol(intensity_data))] %>% set_rownames(., intensity_data$Protein)
}

combine_exprs <- function(exprs_list) {
  Reduce(function(x, y) full_join(x, y, by = "Protein"), exprs_list)
}

clean_expression_data <- function(data) {
  data <- na.omit(data)
  data[data == 0] <- NA
  data <- na.omit(data)
  rownames(data) <- data$Protein
  log2(data[, -1])
}

# ----------------- Step 2: Batch Correction (Ground Truth) -----------------
create_ground_truth <- function(data, batch_factor) {
  sva::ComBat(data, batch = batch_factor)
}

# ------------------ Step 3: Simulate MVs into Ground Truth ------------------
simulate_mvs <- function(data, mv_prop, simbeam = FALSE, batch_factor = NULL, seed=NULL, sample.MNAR = T){
  if(!is.null(seed)){
    set.seed(seed)
  }
  Actual_data = data
  if(simbeam == TRUE){
    batches = unique(batch_factor)
    beams_data = c()
    for(i in batches){
      batch_data = apply(data[,which(batch_factor == i)], 1, mean, na.rm=TRUE)
      beams_data = cbind(beams_data, batch_data)
    }
    colnames(beams_data) = seq(length(batches))
    data = beams_data
  }
  D = as.matrix(data)
  n = length(D)
  mnar_ratio = mv_prop[2]
  mv_prop = mv_prop[1]
  alpha = mv_prop
  beta = mnar_ratio #ratio relative to mv_prop
  if (sample.MNAR == F){
    tmat = matrix(rnorm(n, mean = quantile(D, alpha, na.rm=TRUE), sd = 0.3), ncol=ncol(D), nrow=nrow(D))
  } else if (sample.MNAR == T){
    tmat <- matrix(NA, ncol=ncol(D), nrow=nrow(D))
    for (col in 1:ncol(tmat)){
      tmat[,col] = rnorm(nrow(tmat), mean = quantile(data[,col], alpha, na.rm=TRUE), sd = 0.3)
    }
  }
  pmat = matrix(rbinom(n, 1, beta), ncol=ncol(D), nrow=nrow(D))
  
  
  # inject MNAR
  diff_mat = D < tmat
  diff_mat[diff_mat == FALSE] = 0
  diff_mat[diff_mat == TRUE] = 1
  pmat[pmat == FALSE] = 0
  pmat[pmat == TRUE] = 1
  
  indicator = diff_mat + pmat
  
  msdata = D
  # Inject MNAR
  msdata[indicator == 2] = NA
  
  
  if (simbeam == TRUE){
    msdata = Actual_data
    for (j in batches){
      msdata[which(indicator[,j] == 2),which(batch_factor == j)] = NA
    }
  }
  
  # Determine MCAR/MAR percentage
  total_mcar = n * alpha * (1-beta)
  remaining_observed = which(is.finite(msdata))
  mcar_inds = sample(remaining_observed, round(total_mcar), replace = FALSE)
  
  # Inject MCAR
  if (simbeam == TRUE){
    remaining_observed = which(is.finite(indicator))
    mcar_inds = sample(remaining_observed, round(total_mcar), replace = FALSE)
    indicator[mcar_inds] = 3
    for (j in batches){
      msdata[which(indicator[,j] == 3),which(batch_factor == j)] = NA
    }
  } else{
    msdata[mcar_inds] = NA
  }
  
  return(list(mcar_ind = mcar_inds,
              mnar_ind = indicator,
              msdata = msdata))
}

filter_high_mv_features <- function(data, threshold = 0.6){
  mv_amt = apply(data, 1, function(x) length(which(is.na(x))))
  mv_perc = (mv_amt/ncol(data))
  drop_inds = which(mv_perc > threshold)
  
  if (length(drop_inds) > 0){
    filtered_data = data[-drop_inds,]
  } else {
    filtered_data = data
  }
  return(filtered_data)
}

# ----------------- Step 4 & 5: Subset + Transfer MV Locations -----------------
transfer_mv_pattern <- function(full_data, mv_mask) {
  subset <- full_data[rownames(mv_mask), ]
  subset[is.na(mv_mask)] <- NA
  subset
}

# ------------------ Step 6: Imputation Methods ------------------
impute_methods <- function(data, methods = c("MICE-PMM", "MICE-norm", "KNN-feature", "KNN-sample", "RF"), 
                           ksamp = 2, kfeat = 10, ntrees = 100, sigma = 0.5) {
  out <- list()
  
  if ("MICE-PMM" %in% methods) {
    out$`MICE-PMM` <- mice(janitor::clean_names(data)) %>% complete(.) %>% as.matrix()
  }
  if ("MICE-norm" %in% methods) {
    out$`MICE-norm` <- mice(janitor::clean_names(data), meth = "norm") %>% complete(.) %>% as.matrix()
  }
  if ("KNN-feature" %in% methods) {
    out$`KNN-feature` <- impute.knn(as.matrix(data), k = kfeat)$data
  }
  if ("KNN-feature2" %in% methods) {
    out$`KNN-feature` <- impute.knn(as.matrix(data), k = kfeat, rowmax=0.99, colmax=0.99)$data
  }
  if ("KNN-sample" %in% methods) {
    out$`KNN-sample` <- do.knn(as.matrix(data), k = ksamp)$imputed_data
  }
  if ("RF" %in% methods) {
    out$RF <- t(missForest(t(data), ntree = ntrees)$ximp)
  }
  if ("QRILC" %in% methods) {
    out$QRILC <- impute.QRILC(data, tune.sigma = sigma)[[1]]
  }
  out
}

# ------------------ Step 7: Correct Imputed Batch+ve ------------------
correct_batch_effects <- function(imputed_list, batch_factor) {
  lapply(imputed_list, function(x) sva::ComBat(x, batch = batch_factor))
}

# ------------------ Step 8: Evaluation ------------------
evaluate_rmse <- function(imputed_list, mis_data, true_data, normalizeRMSE = T, batch_rmse = F, batch_factor=NULL) {
  imputed_list <- imputed_list[names(imputed_list) != "True"]
  if(batch_rmse == T){
    rmse_list <- lapply(imputed_list, function(x) Batch_RMSE(x, mis_data, true_data, batch = batch_factor, normalized = normalizeRMSE))
  } else {
    rmse_list <- lapply(imputed_list, function(x) Rmse(x, mis_data, true_data, norm = normalizeRMSE)[1])
  }
  do.call(rbind, rmse_list)
}

evaluate_batch_median <- function(dataset_list, class_factor, batch_factor, set, iter) {
  out <- lapply(dataset_list, function(x) {
    pca <- prcomp(t(x), scale = TRUE, center = TRUE)$x
    batch_median(pca, class_factors = class_factor, batch_factors = batch_factor)
  }) %>% unlist()
  data.frame(value = out, variable = names(out), set = set, iteration = iter)
}

evaluate_dea <- function(dataset_list, class_factor, set, iter) {
  res <- lapply(names(dataset_list), function(name) {
    scores <- limmaDEA(dataset_list[[name]], class_factor, applyfc = TRUE, set.thres = 0.1, output = "DEA")
    data.frame(Method = name, Metric = names(scores), Value = as.numeric(scores))
  })
  do.call(rbind, res) %>%
    mutate(set = set, iteration = iter)
}

evaluate_logfc_error <- function(Set1_dfs, Set2_dfs, Set1_true, class_factor, i) {
  true_logfc <- limmaDEA(Set1_true, class_factor, applyfc = FALSE, set.thres = 0.2, output = "logfc")
  true_logfc <- setNames(as.numeric(true_logfc), rownames(Set1_true)[1:length(true_logfc)])
  DAP <- intersect(names(true_logfc), grep("YEAST|ECOLI", names(true_logfc), value = TRUE))
  
  Set1_logfc <- lapply(Set1_dfs, function(x) {
    out <- limmaDEA(x, class_factor, applyfc = FALSE, set.thres = 0.2, output = "logfc")
    setNames(as.numeric(out), rownames(x)[1:length(out)])
  })
  
  Set2_logfc <- lapply(Set2_dfs, function(x) {
    out <- limmaDEA(x, class_factor, applyfc = FALSE, set.thres = 0.2, output = "logfc")
    setNames(as.numeric(out), rownames(x)[1:length(out)])
  })
  
  err1 <- lapply(Set1_logfc, function(x) abs(x[DAP] - true_logfc[DAP])) %>%
    do.call(rbind, .) %>%
    melt() %>%
    mutate(set = "batch-ve", iteration = i)
  
  err2 <- lapply(Set2_logfc, function(x) abs(x[DAP] - true_logfc[DAP])) %>%
    do.call(rbind, .) %>%
    melt() %>%
    mutate(set = "batch+ve", iteration = i)
  
  rbind(err1, err2)
}

calculate_batch_rmse <- function(true_data, test_data, na.rm = TRUE) {
  # Check if dimensions match
  if (!all(dim(true_data) == dim(test_data))) {
    stop("Datasets must have the same dimensions.")
  }
  
  # Compute squared differences
  sq_diff <- (true_data - test_data)^2
  
  # Compute RMSE
  rmse <- sqrt(mean(sq_diff, na.rm = na.rm))
  
  return(rmse)
}

# ----------------- Run Full BRAID Workflow -----------------
run_braid <- function(df, batch_f, class_f, mv_prop, iterations, mvi_methods = c("MICE-PMM", "MICE-norm", "KNN-feature", "KNN-sample", "RF"),
                      ksamp = 2, kfeat = 10, ntrees = 100, sigma = 0.5) {
  complete_exprs_log <- df
  
  batch_factor <- as.factor(batch_f)
  class_factor <- as.factor(class_f)
  pdat <- data.frame(batch = batch_factor, class = class_factor)
  pdat$FullRunName <- colnames(complete_exprs_log)
  
  start_ground_truth <- create_ground_truth(complete_exprs_log, batch_factor)
  
  save_RMSE <- c()
  save_batchmedian <- c()
  save_DEA <- c()
  save_logfc_error <- c()
  
  for (i in 1:iterations) {
    mv_sim <- simulate_mvs(start_ground_truth, mv_prop = mv_prop, batch_factor = batch_factor, seed = 123456 + i)
    Set1_mis <- filter_high_mv_features(mv_sim$msdata, 0.6)
    Set1_true <- start_ground_truth[rownames(Set1_mis), ]

    Set2_true <- start_ground_truth[rownames(Set1_mis), ]
    Set2_mis <- transfer_mv_pattern(complete_exprs_log, Set1_mis)
    
    Set1_imputed <- impute_methods(Set1_mis, methods = mvi_methods, ksamp = ksamp, kfeat = kfeat, ntrees = ntrees, sigma = sigma)
    Set2_imputed <- impute_methods(Set2_mis, methods = mvi_methods, ksamp = ksamp, kfeat = kfeat, ntrees = ntrees, sigma = sigma)
    Set2_corrected <- correct_batch_effects(Set2_imputed, batch_factor)
    Set2_corrected[["True"]] <- Set2_true
    
    Set1_imputed[["True"]] <- Set1_true
    
    rmse1 <- evaluate_rmse(Set1_imputed, Set1_mis, Set1_true, batch_rmse = F)
    rmse2 <- evaluate_rmse(Set2_corrected, Set2_mis, Set2_corrected$True, batch_rmse = F)
    save_RMSE <- rbind(save_RMSE, cbind(melt(rmse1), set = "batch-ve", iteration = i), cbind(melt(rmse2), set = "batch+ve", iteration = i))
    
    bmed1 <- evaluate_batch_median(Set1_imputed, class_factor, batch_factor, "batch-ve", i)
    bmed2 <- evaluate_batch_median(Set2_corrected, class_factor, batch_factor, "batch+ve", i)
    save_batchmedian <- rbind(save_batchmedian, bmed1, bmed2)
    
    Set1_imputed <- lapply(Set1_imputed, function(x) set_rownames(x, rownames(Set1_true)))
    Set2_corrected <- lapply(Set2_corrected, function(x) set_rownames(x, rownames(Set2_true)))
    
    dea1 <- evaluate_dea(Set1_imputed, class_factor, "batch-ve", i)
    dea2 <- evaluate_dea(Set2_corrected, class_factor, "batch+ve", i)
    save_DEA <- rbind(save_DEA, dea1, dea2)
    
    logfc_error <- evaluate_logfc_error(Set1_imputed, Set2_corrected, Set1_true, class_factor, i)
    save_logfc_error <- rbind(save_logfc_error, logfc_error)
    
  }
  
  list(RMSE = save_RMSE, batchmedian = save_batchmedian, DEA = save_DEA, logfc_error = save_logfc_error)
}


# --------------------------------------------------
# ----------------- BRAID example usage: -----------------
# --------------------------------------------------
vp_wd <- paste0(wd,"/BRAID_BEGONE/Van Puyvelde dataset/")
design_files <- c("HYE5600735_LFQ_FragPipe_design.tsv", "HYE6600735_LFQ_FragPipe_design.tsv", "HYEqe735_LFQ_FragPipe_design.tsv", "HYEtims735_LFQ_FragPipe_design.tsv")
intensity_files <- c("HYE5600735_LFQ_FragPipe_pro_intensity.tsv", "HYE6600735_LFQ_FragPipe_pro_intensity.tsv", "HYEqe735_LFQ_FragPipe_pro_intensity.tsv", "HYEtims735_LFQ_FragPipe_pro_intensity.tsv")


# Prepare the data
design_list <- load_design_data(vp_wd, design_files)
intensity_list <- load_intensity_data(vp_wd, intensity_files)
exprs_list <- lapply(intensity_list, extract_exprs)
combined_exprs <- combine_exprs(exprs_list)
input_df <- clean_expression_data(combined_exprs)

batch_f <- as.factor(rep(1:4, times = sapply(design_list, nrow)))
class_f <- do.call(rbind, design_list) %>% .$condition %>% as.factor() %>% as.numeric() %>% as.factor()

results <- run_braid(input_df, batch_f = batch_f, class_f = class_f, mv_prop = c(0.3, 0.7), iterations = 5, mvi_methods = c("KNN-feature", "KNN-sample"), ksamp = 4)

# Visualization can be added here using ggplot2 based on `results`
plot_RMSE <- results$RMSE
plot_RMSE$Var1 = as.factor(plot_RMSE$Var1)
plot_RMSE$value = as.numeric(plot_RMSE$value)
ggplot(plot_RMSE, aes(x=Var1, y=value, fill=set)) + geom_boxplot(position="dodge", outliers = T) +
  labs(x="", y="NRMSE", fill="") +
  theme(axis.text.x = element_text(angle=20, hjust=1, vjust=1))

plot_logfc_error <- results$logfc_error
plot_logfc_error$Var1 = as.factor(plot_logfc_error$Var1)
plot_logfc_error$value = as.numeric(plot_logfc_error$value)
ggplot(plot_logfc_error, aes(x=Var1, y=value, fill=set)) + 
  geom_boxplot(position="dodge", outliers=T) +
  labs(x="MVI method", y="Absolute LogFC error", fill="") +
  theme(axis.text.x = element_text(angle=20, hjust=1, vjust=1))


# ========== BEGONE: MVI Robustness Evaluation ==========
# ------------------ Step 1: Generate ground-truth ------------------
# Perform BRAID framework up to step 3

# ------------------ Step 2: Simulate Progressive Batch Effects ------------------
sim_batch_effects <- function(data, 
                              batch_assignments, 
                              class_assignments = NULL,
                              reference_batch = NULL,
                              shift = 0.5, scale = 0.1,
                              add.sd = shift/5, multi.sd = scale/5,
                              additive = T, multiplicative = T, 
                              class_heterogeneous = F, # Simulates heterogeneous batch effects for different classes in the same batch. Requires class assignments.
                              seed = 1234567890, parameters = F) {
  set.seed(seed)
  # Check that batch_assignments match the number of columns
  if (length(batch_assignments) != ncol(data)) {
    stop("Length of batch_assignments must match the number of columns in the data.")
  }
  if (abs(scale) > 1.5){
    stop("Scale value cannot exceed 1.5.")
  }
  
  if (is.null(colnames(data))){
    colnames(data) <- seq(1:ncol(data))
  }
  
  if (!is.null(reference_batch)){
    if (!(reference_batch %in% unique(batch_assignments))){
      stop("Reference batch not found within batch assignments.")
    }
    # Original assignments
    store_batch_assignments = batch_assignments
    store_class_assignments = class_assignments
    
    # Store reference batch data
    ref_batch = data[,batch_assignments == reference_batch]
    
    # New batch and class assignments
    batch_assignments = batch_assignments[batch_assignments != reference_batch]
    class_assignments = class_assignments[batch_assignments != reference_batch]
  }
  # Unique batch labels
  batches <- unique(batch_assignments)
  
  # Duplicate input data to introduce batch effects
  data_with_effects <- data
  if (!is.null(reference_batch)){
    data_with_effects <- data[,store_batch_assignments != reference_batch]
  }
  # Apply batch effects
  loc_shift <- c()
  scale_shift <- c()

  loc_vecs <- list()
  scale_vecs <- list()
  
  # Turn off batch effect introduction if parameters are 0
  if (shift == 0){
    additive = F
  }
  if (scale == 0){
    multiplicative = F
  }
  if (additive == F & multiplicative == F){
    if (parameters == T){
      cat("No parameters are available as no batch effects were added")
      if (!is.null(reference_batch)){
        return(list(batch_mat = data,
                    location_shift = NULL,
                    location_effects = matrix(0, nrow=nrow(data), ncol=length(batches)),
                    scale_shift = NULL,
                    scale_effects = matrix(1, nrow=nrow(data), ncol=ncol(data))
                    # noise = epsilon
        )
        )
      }
      return(list(batch_mat = data_with_effects,
                  location_shift = loc_shift,
                  location_effects = matrix(0, nrow=nrow(data), ncol=length(batches)),
                  scale_shift = scale_shift,
                  scale_effects = matrix(1, nrow=nrow(data), ncol=ncol(data))
                  # noise = epsilon
      )
      )
    } else {
      if (!is.null(reference_batch)){
        return(data)
      }
      return(data_with_effects)
    }
  }
  for (batch in batches) {
    # Get indices for the current batch
    batch_indices <- which(batch_assignments == batch)
    
    # Generate batch effect vector
    location_effect_vec <- 0
    scale_effect_vec <- 1
    
    if (additive == T){
      location_effect <- runif(1, min = -(shift/2), max = (shift/2))
      if (class_heterogeneous == F){
        location_effect_vec <- rnorm(nrow(data), mean = location_effect, sd = add.sd)
        loc_shift <- c(loc_shift, location_effect) # save mean shift
        loc_vecs[[batch]] <- location_effect_vec # save mean shift effect vector
      } else if (class_heterogeneous == T){
        # Create a list with length(unique(classes)) to store batch effect vectors
        class_location_effect_list <- as.list(unique(class_assignments))
        for (class in unique(class_assignments)){
          location_effect_vec <- rnorm(nrow(data), mean = location_effect, sd = add.sd)
          class_location_effect_list[[class]] <- location_effect_vec
        }
        loc_shift <- c(loc_shift, location_effect) # save mean shift
        loc_vecs[[batch]] <- class_location_effect_list # save mean shift effect vector
      }
    }
    if (multiplicative == T){
      scale_effect <- runif(1, min = 1-scale, max = 1+scale)
      if (class_heterogeneous == F){
        scale_effect_vec <- rnorm(nrow(data), scale_effect, sd=multi.sd)
        scale_shift <- c(scale_shift, scale_effect)
        scale_vecs[[batch]] <- scale_effect_vec
      } else if (class_heterogeneous == T){
        # Create a list with length(unique(classes)) to store batch effect vectors
        class_scale_effect_list <- as.list(unique(class_assignments))
        for (class in unique(class_assignments)){
          scale_effect_vec <- rnorm(nrow(data), scale_effect, sd=multi.sd)
          class_scale_effect_list[[class]] <- scale_effect_vec
        }
        scale_shift <- c(scale_shift, scale_effect) # save scale shift
        scale_vecs[[batch]] <- class_scale_effect_list # save scale shift effect vector
      }
    }
    # Apply effects: Shift and scale
    if (class_heterogeneous == F){
      data_with_effects[, batch_indices] <- (data_with_effects[, batch_indices] + location_effect_vec) * scale_effect_vec
    } else if (class_heterogeneous == T){
      for (class in class_assignments){
        class_indices = which(class_assignments == class)
        batch_class_indices = class_indices[which(class_indices %in% batch_indices)]
        if (additive == T & multiplicative == T){
          data_with_effects[, batch_class_indices] <- (data_with_effects[, batch_class_indices] + class_location_effect_list[[class]]) * class_scale_effect_list[[class]]
        } else if (additive == F & multiplicative == T){
          data_with_effects[, batch_class_indices] <- data_with_effects[, batch_class_indices] * class_scale_effect_list[[class]]
        } else if (additive == T & multiplicative == F){
          data_with_effects[, batch_class_indices] <- data_with_effects[, batch_class_indices] + class_location_effect_list[[class]]
        } 
      }
    }

  }
  
  if (!is.null(reference_batch)){
    data_with_effects = cbind(data_with_effects, ref_batch)
    data_with_effects = data_with_effects[,colnames(data)]
  }
  
  if (parameters == T){
    if (class_heterogeneous == F){
      return(list(batch_mat = data_with_effects,
                  location_shift = loc_shift,
                  location_effects = as.data.frame(loc_vecs),
                  scale_shift = scale_shift,
                  scale_effects = as.data.frame(scale_vecs)
                  # noise = epsilon
      )
      )
    } else if (class_heterogeneous == T){
      return(list(batch_mat = data_with_effects,
                  location_shift = loc_shift,
                  location_effects = loc_vecs,
                  scale_shift = scale_shift,
                  scale_effects = scale_vecs
                  # noise = epsilon
      )
      )
    }
  }
  return(data_with_effects)
}

# ------------------ Step 3: Imputation and batch correction ------------------
# Use steps 6 & 7 from BRAID

# ------------------ Step 4: Evaluation ------------------
# Use step 8 from BRAID

# ------------------ Step 5: Run several iterations ------------------
#' Note:
#' This function is only an example to run BEGONE for your dataset.
#' You should customize the function to suit your needs.
run_begone <- function(df, batch_f, class_f, mv_prop,
                       n_iter = 10,
                       mean_vec = seq(0, by = 0.5, length.out = 10),
                       sd_vec = seq(0, by = 0.05, length.out = 10),
                       additive_sd = 0,
                       multiplicative_sd = 0,
                       sample_features = F, feature_subsample = 0.1,
                       simbeams = F,
                       mvi_methods = c("MICE-PMM", "MICE-norm", "KNN-feature", "KNN-sample", "RF"),
                       seed = 1234567890)
{
  set.seed(seed)
  
  complete_exprs_log <- df
  
  batch_factor <- as.factor(batch_f)
  class_factor <- as.factor(class_f)
  pdat <- data.frame(batch = batch_factor, class = class_factor)
  
  start_ground_truth <- create_ground_truth(complete_exprs_log, batch_factor)
  
  feature_count <- nrow(df)
  all_datasets <- list()
  
  for (mean in 1:length(mean_vec)){
    for (sd in 1:length(sd_vec)){
      for (i in 1:n_iter) {
        cat("Running iteration:", i, "\n")
        cat("Running iteration", i, "for:", "Mean", mean_vec[mean], ", SD", sd_vec[sd])
        if (sample_features == T){
          sampled_features <- sample(1:feature_count, size = ceiling(feature_count * feature_subsample))
          df <- df[, sampled_features]
        }
        # Simulate MVs and filter high MV features
        mv_sim <- simulate_mvs(start_ground_truth, mv_prop, batch_factor, simbeam = simbeams, seed + i)
        Set1_mis <- filter_high_mv_features(mv_sim$msdata, 0.6)
        Set1_true <- start_ground_truth[rownames(Set1_mis), ]
        
        # Add batch effects
        batch_affected <- sim_batch_effects(Set1_mis, batch_assignments=batch_factor, class_assignments=class_factor, shift = mean_vec[mean], scale = sd_vec[sd], add.sd=additive_sd, multi.sd=multiplicative_sd, class_heterogeneous = F, reference_batch = unique(batch_f)[1])
        
        imputed <- impute_methods(batch_affected, methods = mvi_methods)
        corrected <- correct_batch_effects(imputed, batch_factor)
        corrected[["True"]] <- Set1_true
        corrected[["Missing"]] <- Set1_mis
        
        mean_label <- paste0("mean_", mean_vec[mean])
        sd_label <- paste0("sd_", sd_vec[sd])
        iter_label <- paste0("iter_", i)
        
        if (is.null(all_datasets[[mean_label]])) all_datasets[[mean_label]] <- list()
        if (is.null(all_datasets[[mean_label]][[sd_label]])) all_datasets[[mean_label]][[sd_label]] <- list()
        
        all_datasets[[mean_label]][[sd_label]][[iter_label]] <- corrected
        
      }
    }
  }
  return(all_datasets)
}


# --------------------------------------------------
# ------------- BEGONE Example usage: --------------
# --------------------------------------------------
# Testing on the Poulos dataset
# Load and prepare the dataset
testdatanew <- read.csv(paste0(wd,"/BRAID_BEGONE/Poulos dataset/Peptide_intensity_matrix_b9369842-f3cf-4383-9e5c-6734dadcfbc9.csv"))
testmetanew <- readxl::read_xlsx(paste0(wd,"/BRAID_BEGONE/Poulos dataset/Mapping_file_PXD015912.xlsx"))
testmetanew <- testmetanew[testmetanew$Filetype == "SWATH",]
testmetanew$Filename <- tolower(testmetanew$Filename)
testmetanew$Filename <- sub("-", ".", testmetanew$Filename)
testmetanew$Filename <- paste0("X",testmetanew$Filename)

# Match columns
match_columns <- intersect(colnames(testdatanew), testmetanew$Filename)
testmetanew <- testmetanew[which(testmetanew$Filename %in% match_columns),]
existing_columns <- which(colnames(testdatanew) %in% testmetanew$Filename)

exprsnew <- testdatanew[,existing_columns]
exprsnew <- cbind(exprsnew, Protein=testdatanew$Protein,
                  Peptide=testdatanew$Peptide)
rownames(exprsnew) <- exprsnew$Peptide
newcols <- strsplit(testmetanew$Information, "_")
newcols <- do.call(rbind, newcols)
colnames(newcols) <- c("Day", "Batch", "Condition")
testmetanew <- cbind(testmetanew, newcols)


# BEGONE functions made specifically for the Poulos dataset
run_begone_poulos <- function(df, testmetanew, mv_prop,
                              mean_vec = seq(0, by = 0.5, length.out = 10),
                              sd_vec = seq(0, by = 0.05, length.out = 10),
                              additive_sd = NULL,
                              multiplicative_sd = NULL,
                              sample_features = FALSE, feature_subsample = 0.1,
                              simbeams = FALSE,
                              mvi_methods = c("MICE-PMM", "MICE-norm", "KNN-feature", "KNN-sample", "RF")) {
  set.seed(1234567890)
  seeds=c()
  all_datasets <- list()
  
  for (mean_val in mean_vec) {
    mean_key <- paste0("mean_", mean_val)
    all_datasets[[mean_key]] <- list()
    
    for (sd_val in sd_vec) {
      # if (sd_val == 0.02){
      #   browser()
      # }
      
      sd_key <- paste0("sd_", sd_val)
      iteration_results <- list()
      
      for (i in 1:length(unique(testmetanew$Day))) {
        
        cat("Running iteration", i, "for Mean =", mean_val, ", SD =", sd_val, "\n")
        
        # Extract 1 day only
        metaday1 <- testmetanew[which(testmetanew$Day == unique(testmetanew$Day)[i]),]
        exprsday1 <- df[,metaday1$Filename]
        
        
        fdata <- df[match(rownames(exprsday1), df$Peptide) ,c("Protein", "Peptide")]
        colnames(fdata) <- c("protein","peptide")
        exprsday1_prot <- aggregate_protein_expression(exprsday1, fdata, is.log=F)
        
        mean(is.na(exprsday1_prot))
        
        # Remove any samples with high MV proportions
        sample_mvs <- apply(exprsday1_prot, 2, function(x) mean(is.na(x)))
        drop_samples <- which(sample_mvs > 0.5)
        if (length(drop_samples) > 0){
          exprsday1 <- exprsday1[,-drop_samples]
          exprsday1_prot <- exprsday1_prot[,-drop_samples]
          metaday1 <- metaday1[-drop_samples,]
        }
        completenew <- na.omit(exprsday1_prot)
        colnames(completenew) <- colnames(exprsday1)
        
        # Extract only samples 1 and 5
        meta_test <- metaday1[which(grepl("Sample 1", metaday1$Condition) | grepl("Sample 5", metaday1$Condition)),]
        exprs_test <- completenew[,meta_test$Filename]
        
        # Set meta data
        batch_factor <- as.factor(meta_test$Batch)
        class_factor <- as.factor(meta_test$Condition)
        pdat <- data.frame(batch=as.factor(meta_test$Batch), class=as.factor(meta_test$Condition))
        batches <- unique(batch_factor)
        
        # if(sd_val == 0.05){
        #   browser()
        # }
        
        
        start_ground_truth <- sva::ComBat(exprs_test, batch_factor)
        df_iter <- start_ground_truth
        
        
        if (sample_features) {
          feature_count <- nrow(df_iter)
          sampled_features <- sample(1:feature_count, size = ceiling(feature_count * feature_subsample))
          df_iter <- df_iter[, sampled_features]
        }
        
        if (is.null(additive_sd)) {
          additive_sd <- mean_val / 5
        }
        if (is.null(multiplicative_sd)) {
          multiplicative_sd <- sd_val / 5
        }
        
        # We use the same seed for each batch effect simulation
        # as the dataset is different each time
        batch_affected <- sim_batch_effects(
          df_iter,
          batch_assignments = batch_factor,
          class_assignments = class_factor,
          shift = mean_val,
          scale = sd_val,
          class_heterogeneous = FALSE,
          reference_batch = "M01",
          parameters = F
        )
        
        start_Set1_true <- sva::ComBat(batch_affected, batch_factor, ref.batch = "M01")
        
        mv_sim <- simulate_mvs(start_Set1_true, mv_prop = c(0.3, 0.7), batch_factor, simbeam = simbeams, seed = 1234567890 + i)
        Set1_mis <- filter_high_mv_features(mv_sim$msdata, 0.6)
        Set1_true <- batch_affected[rownames(Set1_mis), ]
        
        Set2_mis <- Set1_true
        Set2_mis[is.na(Set1_mis)] <- NA
        
      
        Set1_true <- sva::ComBat(Set1_true, batch_factor, ref.batch = "M01")

        imputed <- impute_methods(Set2_mis, methods = mvi_methods)
        corrected <- lapply(imputed, function(x) sva::ComBat(x, batch_factor, ref.batch = "M01"))
        names(imputed) <- paste0("Uncorrected_", names(imputed))
        corrected <- c(imputed, corrected)
        corrected[["True"]] <- Set1_true
        corrected[["Missing"]] <- Set2_mis
        
        iteration_results[[paste0("iter_", i)]] <- corrected
        
        seeds <- rbind(seeds, (.Random.seed[1:10]))
      }
      all_datasets[[mean_key]][[sd_key]] <- iteration_results
    }
  }
  
  return(list(res=all_datasets, seeds=seeds))
}

evaluate_logfc_error_poulos <- function(imputed, true, class_factor, i) {
  true_logfc <- limmaDEA(true, class_factor, output = "logfc")
  true_logfc <- setNames(as.numeric(true_logfc), rownames(true)[1:length(true_logfc)])
  DAP <- intersect(names(true_logfc), grep("YEAST|ECOLI", names(true_logfc), value = TRUE))
  
  imputed_logfc <- limmaDEA(imputed, class_factor, output = "logfc")
  imputed_logfc <- setNames(as.numeric(imputed_logfc), rownames(imputed)[1:length(imputed_logfc)])
  
  err1 <- abs(imputed_logfc[DAP] - true_logfc[DAP])
  
  err1
}


prepare_begone_results_poulos <- function(all_datasets, meta) {
  results_list <- list()
  
  for (mean_key in names(all_datasets)) {
    for (sd_key in names(all_datasets[[mean_key]])) {
      iter_block <- all_datasets[[mean_key]][[sd_key]]
      method_scores <- list()
      method_error <- list()
      batch_scores <- list()
      rowvars <- list()
      for (iter_idx in names(iter_block)) {
        metaday1 <- meta[which(meta$Day == unique(meta$Day)[as.numeric(strsplit(iter_idx, "_")[[1]][2])]),]
        meta_test <- metaday1[which(grepl("Sample 1", metaday1$Condition) | grepl("Sample 5", metaday1$Condition)),]
        
        
        
        iter_data <- iter_block[[iter_idx]]
        true_data <- iter_data[["True"]]
        # missing_mask <- is.na(iter_data[["Missing"]])
        missing_mask <- iter_data[["Missing"]]
        
        
        meta_test <- meta[meta$Filename %in% colnames(true_data),]
        
        class_f <- meta$Condition[match(colnames(true_data), meta$Filename)]
        
        for (method in setdiff(names(iter_data), c("True", "Missing"))) {
          imputed <- iter_data[[method]]
          nrmse <- Rmse(imputed, missing_mask, true_data, norm=T)[[1]]
          method_scores[[method]] <- c(method_scores[[method]], nrmse)
          
          logfc_err <- evaluate_logfc_error_poulos(imputed, iter_data[["True"]], class_factor = class_f)
          method_error[[method]] <- c(method_error[[method]], logfc_err)
          # browser()
          batch_rmse <- mean(Batch_RMSE(imputed, missing_mask, true_data, batch = meta_test$Batch))
          batch_scores[[method]] <- c(batch_scores[[method]], batch_rmse)
          
          rowvars[[method]] <- c(rowvars[[method]], apply(imputed, 1, var, na.rm=T))
        }
      }
      
      for (method in names(method_scores)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          method = method,
          mean = as.numeric(gsub("mean_", "", mean_key)),
          sd = as.numeric(gsub("sd_", "", sd_key)),
          NRMSE = mean(method_scores[[method]]),
          NRMSE_sd = sd(method_scores[[method]]),
          LogFC_error = median(method_error[[method]]),
          batch_rmse = mean(batch_scores[[method]]),
          feature_var = mean(rowvars[[method]])
        )
      }
    }
  }
  
  results_df <- do.call(rbind, results_list)
  return(results_df)
}


poulos_time_s <- Sys.time()
test_begone_poulos <- run_begone_poulos(exprsnew, testmetanew = testmetanew, mv_prop=c(0.3, 0.7),
                                                mean_vec = seq(0, by = 0.5, length.out = 10),
                                                sd_vec = seq(0, by = 0.05, length.out = 10),
                                                feature_subsample = F,
                                                mvi_methods = c("KNN-feature", "KNN-sample")
)
poulos_time_e <- Sys.time()
poulos_time <- poulos_time_e - poulos_time_s
print(poulos_time)

begone_results_poulos <- prepare_begone_results_poulos(test_begone_poulos$res, meta=testmetanew)

res <- begone_results_poulos

# Plotting the results
begone_nrmse_plots <- list()

for(method in unique(res$method)){
  begone_nrmse_plots[[method]] <- ggplot(res[res$method==method,], aes(x=mean, y=sd, fill=NRMSE)) + 
    geom_tile() +
    labs(title=method, fill="NRMSE", x="Additive", y="Multiplicative") +
    scale_fill_gradient(low = "#8Fd9FB", high = "#022333") +
    theme_minimal() +
    theme(panel.grid=element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.text.y = element_text(face="bold"),
          legend.text = element_text(face="bold"),
          title = element_text(face="bold"))
}
ggarrange(begone_nrmse_plots[[3]],
          begone_nrmse_plots[[4]])
