# KNN-sample function
do.knn<-function(df, k, weighted =TRUE, distance="euclidean", percent=0.75){
  # imp.knn <- df
  # imp.knn[is.finite(df) == FALSE] <- NA
  mv_ind<-which(is.na(df), arr.ind = TRUE)
  mv_rows<-unique(mv_ind[,1])
  
  klist<-list()
  
  if(distance == "euclidean"){
    dist_mat = as.matrix(dist(t(df), method="euclidean"))
    imp.knn <- df
    imp.knn[is.finite(df) == FALSE] <- NA
  } else if(distance == "truncation"){
    ParamMat <- EstimatesComputation(df, perc = percent)
    
    means<-ParamMat[,1]
    sd<-ParamMat[,2]
    df<-(df-means)/sd
    
    imp.knn <- df
    imp.knn[is.finite(df) == FALSE] <- NA
    
    dist_mat = as.matrix(cor(df, use = "pairwise.complete.obs"))
  }
  for (i in mv_rows){
    klist[[i]]<-list()
    feature_mvs<-which(is.na(df[i,]))
    feature_obs<-which(!is.na(df[i,]))
    # cand_vectors<-df[,feature_obs]
    for (j in 1:length(feature_mvs)){
      feature_mvs_ind<-feature_mvs[j]
      dist <- dist_mat[feature_mvs_ind, feature_obs] # Get distances for every observed sample in feature
      
      if(distance == "truncation"){
        r <- dist
        rdist <- 1 - abs(dist)
        rdist[is.nan(rdist) | is.na(rdist)] <- Inf
        rdist[rdist==0]<-ifelse(is.finite(min(rdist[rdist>0])), min(rdist[rdist>0])/2, 1)
        rdist[abs(r) == 1] <- Inf
        dist <- rdist
      }
      
      # #if (length(feature_obs) < k){
      # #  imp.knn[i,feature_mvs_ind]<-mean(df[,feature_mvs_ind], na.rm = TRUE)
      # #} else {
      # calc.dist<-data.frame((cand_vectors-target_vector)^2)
      # dist<-sqrt(colMeans(calc.dist, na.rm = TRUE))
      # dist[is.nan(dist) | is.na(dist)] <- Inf
      # dist[dist==0] <- ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
      
      if (length(dist) < k) {
        stop(message = "Fewer than K finite distances found")
        #imp.knn[i,feature_mvs_ind]<-mean(df[i,], na.rm = TRUE)
      } else {
        k_sample_ind <- order(dist)[1:k]
        k_samples <- feature_obs[k_sample_ind]
        if(distance == "euclidean"){
          wghts <- 1/dist[k_sample_ind]/sum(1/dist[k_sample_ind])
        } else if(distance == "truncation"){
          wghts <- (1/dist[k_sample_ind]/sum(1/dist[k_sample_ind])) * sign(r[k_sample_ind])
        }
        if (weighted == TRUE){
          imp.knn[i, feature_mvs_ind] <- wghts %*% df[i, k_samples]
        } else if (weighted == FALSE){
          imp.knn[i, feature_mvs_ind] <- mean(df[i, k_samples])
        }
      }
      #}
      klist[[i]][[feature_mvs_ind]]<-k_samples
    }
  }
  
  if(distance=="truncation") {
    imp.knn <- (imp.knn * sd) + means
  }
  return(list(original_data=df, imputed_data=imp.knn, klist=klist))
  #return(imp.knn)
}

# Limma differential expression analysis wrapper function
limma.wrapper <- function(df, class, fcthres=T, thres=0.5){
  design = model.matrix(~as.factor(class))
  colnames(design) <- c("Intercept", "Contrast")
  lmtest = lmFit(df, design=design)
  lmtest <- eBayes(lmtest)
  pvals = lmtest$p.value[,"Contrast"]
  qvals = p.adjust(pvals, method="fdr")
  if (fcthres == T){
    lmtest2 = treat(lmtest, lfc=thres)
    topTreat(lmtest2,coef=2)
    pvals = lmtest2$p.value[,"Contrast"]
    qvals = p.adjust(pvals, method="fdr")
  } else {
    lmtest2 <- lmtest
  }
  lmtest2$q.value <- qvals
  return(lmtest2)
}

limmaDEA <- function(df, class_factor, applyfc=T, set.thres=0.5, output="features"){
  # browser()
  
  ## perform limma
  a=limma.wrapper(df, class=as.factor(class_factor), fcthres=applyfc, thres=set.thres)
  
  # get logfc estimate
  logfc <- a$coefficients[,2]
  
  result=data.frame(pval = a$p.value[,2], qval = a$q.value, logfc = logfc)
  
  rownames(result)=rownames(df)
  result = cbind(rownames(df), result)
  colnames(result)=c("protein", "pval", "qval", "logfc")
  
  # save result to list
  DA_protein_inds = which(result[,"qval"]<0.05)
  
  if (output == "all"){
    out <- result
  }
  if (output == "significant"){
    out <- result$protein[DA_protein_inds]
  }
  if (output == "pval"){
    out <- result[,'pval']
  }
  if (output == "qval"){
    out <- result[,'qval']
  }
  if (output == "logfc"){
    out <- result[,'logfc']
  }
  
  if (output == "roc"){
    # Get ROC data
    rocdata = as.data.frame(result)
    rocdata$trueDA = FALSE
    rocdata$trueDA[which(grepl('YEAST',rocdata$protein))] <- TRUE
    rocdata$trueDA[which(grepl('ECOLI',rocdata$protein))] <- TRUE
    
    # sort by increasing q-value
    rocdata = rocdata[order(rocdata$qval),]
    
    truehits<-c()
    count=0
    for (j in 1:length(rocdata$trueDA)){
      if (rocdata$trueDA[j]==TRUE){
        count=count+1
      }
      truehits=c(truehits,count)
    }
    
    TPR<-c()
    Precision<-c()
    FPcount<-0
    for (h in 1:length(truehits)){
      # TPR
      val<-(truehits[h])/(max(truehits))
      TPR<-c(TPR,val)
      
      # Precision
      if (rocdata$trueDA[h]==FALSE){
        FPcount<-FPcount+1
      }
      Precision<-c(Precision, (truehits[h]/(truehits[h]+FPcount)))
    }
    
    # FPR
    FPR<-c()
    FPcount<-0
    for (g in 1:length(rocdata$trueDA)){
      if (rocdata$trueDA[g]==FALSE){
        FPcount<-FPcount+1
      }
      FPR<-c(FPR,(FPcount/(length(which(rocdata$trueDA==FALSE)))))
    }
    
    response<-rocdata$trueDA
    tups<-which(rocdata$trueDA == 'TRUE')
    response[tups]<-1
    response[!tups]<-0
    
    rocdata <- cbind(rocdata,truehits=truehits,TPR=TPR,FPR=FPR,Precision=Precision,response=response)
    out = rocdata
  }
  if(output == "DEA"){
    # generate confusion_matrix
    Positives = result$protein[DA_protein_inds]
    TP <- length(c(Positives[which(grepl("YEAST", Positives))],
                   Positives[which(grepl("ECOLI", Positives))]))
    FP <- length(Positives[which(grepl("HUMAN", Positives))])
    Negatives = result$protein[-DA_protein_inds]
    TN <- length(Negatives[which(grepl("HUMAN", Negatives))])
    FN <- length(c(Negatives[which(grepl("YEAST", Negatives))],
                   Negatives[which(grepl("ECOLI", Negatives))]))
    
    scores = c(TP=TP,FP=FP,TN=TN,FN=FN)
    if (any(!is.finite(scores))){
      scores[which(!is.finite(scores))] = 0
      TP = scores["TP"]
      FP = scores["FP"]
      TN = scores["TN"]
      FN = scores["FN"]
    }
    TPR = TP/(TP + FN)
    FPR = FP/(FP + TN)
    Precision = TP/(TP + FP)
    Fscore = (2 * Precision * TPR)/(Precision + TPR)
    
    out <- c(TPR = TPR,
             FPR = FPR,
             Precision = Precision,
             Fscore = Fscore)
  }
  
  # get output
  return(out)
}

# Calculate RMSE
Rmse <- function(imp, mis, true, norm = FALSE){
  imp <- as.matrix(imp)
  mis <- as.matrix(mis)
  true <- as.matrix(true)
  missIndex <- which(is.na(mis))
  errvec <- imp[missIndex] - true[missIndex]
  rmse <- sqrt(mean(errvec^2))
  if (norm) {
    rmse <- rmse/sd(true[missIndex])
  }
  return(list(rmse,sd(true[missIndex])))
}

# Within batch RMSE function
Batch_RMSE <- function(imp, mis, true, batch, normalized = T){
  batch_factors <- as.factor(batch)
  unique_batches <- unique(batch_factors)
  output <- c()
  for(i in unique_batches){
    temp_imp <- imp[,batch_factors == i]
    temp_mis <- mis[,batch_factors == i]
    temp_true <- true[,batch_factors == i]
    RMSE_res <- Rmse(temp_imp, temp_mis, temp_true, norm = normalized)[[1]]
    output <- c(output, RMSE_res)
  }
  return(output)
}

# PCA-based batch effect quantification
batch_median <- function(df, class_factors, batch_factors){
  # pcadf = PCA matrix, col = PCs, rows = samples
  # class_factors = vector (?) of class factors
  # batch_factors = vector (?) of batch factors
  # adapted from Limsoon's method
  # for each class, between possible batch combinations of 2s,
  # the squared distance is calculated and summed per PC
  # i.e. for pc in PCn, for c in n_classes, and for i,j in batches where i != j, and i > j
  # batch_median(pc) = (median_pc,c,i - median_pc,c,j)
  
  # get unique class and batch factors
  batch_factors = as.numeric(as.factor(batch_factors))
  class_factors = as.numeric(as.factor(class_factors))
  unique_batch = unique(batch_factors)
  unique_class = unique(class_factors)
  #output vector
  out_median_sum = c()
  
  # for PC in PCn
  for (pc in 1:dim(df)[1]){
    # distance per PC
    distance = 0.0
    # for combinat package, when vector length == 2 and m == 2, it will return 1x2 vector.
    # so a separate code block is used made if unique batches == 2
    if (length(unique_batch) == 2){
      # for class in unique classes
      for (cl in unique_class){
        # b1_indices/b2_indices are indices of the batch_factors vector for batch 1 and batch 2 respectively
        b1_indices = intersect(which(batch_factors %in% c(unique_batch[1])),which(class_factors %in% c(cl)))
        b2_indices = intersect(which(batch_factors %in% c(unique_batch[2])),which(class_factors %in% c(cl)))
        # b1_val and b2_val are medians of b1 and b2 respectively
        b1_val = median(df[b1_indices, pc])
        b2_val = median(df[b2_indices, pc])
        # cl_dist is the final (x1 - x2)**2 distance
        cl_dist = (b1_val-b2_val)**2
        # since the values are by PC, the distances for all classes are summed per PC
        distance = distance + cl_dist
        
      }
    }
    # else clause refers to if n batch factors > 2
    else{
      # combinations for x batch factors with no replacement for 2 permutations
      batch_combi = combn(1:length(unique_batch),2)
      # distnace per PC
      distance = 0.0
      # for each batch permutations & for each class
      for (perm in 1:dim(batch_combi)[2]){
        for (cl in unique_class){
          # b1/b2 = batch factor 1 and batch factor 2
          b1 = batch_combi[,perm][1]
          b2 = batch_combi[,perm][2]
          # # b1_indices/b2_indices are indices of the batch_factors vector for batch 1 and batch 2 respectively
          b1_indices = intersect(which(batch_factors %in% c(b1)),which(class_factors %in% c(cl)))
          b2_indices = intersect(which(batch_factors %in% c(b2)),which(class_factors %in% c(cl)))
          # b1_val and b2_val are medians of b1 and b2 respectively
          b1_val = median(df[b1_indices, pc])
          b2_val = median(df[b2_indices, pc])
          # cl_dist is the final (x1 - x2)**2 distance
          cl_dist = (b1_val-b2_val)**2
          # since the values are by PC, the distances for all classes are summed per PC
          distance = distance + cl_dist
        }
      }
    }
    
    #at the end of each PC loop, each sum is appended, where n distances == n PCs 
    out_median_sum = c(out_median_sum,distance)
  }
  return(sum(out_median_sum))
}

# Aggregate peptide level expression to protein level expression
aggregate_protein_expression <- function(peptide_data, feature_data, is.log=TRUE) {
  # Ensure row names in feature_data match column names in peptide_data
  if (!all(rownames(peptide_data) %in% feature_data$peptide)) {
    stop("Peptide IDs in feature_data do not match column names in peptide_data.")
  }
  peptide_data <- apply(peptide_data, 2, as.numeric)
  
  if (is.log == TRUE){
    peptide_data = 2^peptide_data
  }
  
  # Merge peptide expression data with feature data to get protein-level expression
  protein_data <- matrix(nrow=length(unique(feature_data[,"protein"])), ncol=0)
  for (i in 1:ncol(peptide_data)){
    helper <- cbind(feature_data, intensity=peptide_data[,i])
    pep_to_prot <- helper %>%
      group_by(protein) %>%
      summarise(total_intensity = sum(intensity, na.rm=TRUE))
    protein_data <- cbind(protein_data, pep_to_prot$total_intensity)
  }
  rownames(protein_data) <- pep_to_prot$protein
  protein_data = log2(protein_data)
  protein_data[protein_data==0] <- NA
  protein_data[is.infinite(protein_data)] <- NA
  return(protein_data)
}

# Function used to organize BEGONE results obtained specifically from the poulos dataset
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
