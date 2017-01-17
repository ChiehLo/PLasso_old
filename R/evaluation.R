evaluation <- function(sample, result){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample);

  RMSE <- RMSE_l0(sample$cor, result$cor);
  ICOV <- ICOV_l0(ginv(sample$cor), result$icov);
  RMSE_f <- RMSE_F(sample$cor, result$cor);
  ICOV_f <- ICOV_F(ginv(sample$cor), result$icov); 
  # A <- sample$adj;
  # A_opt <- result$adj;
  # #A_opt <- result$icov;
  # perf <- ROC(A, A_opt, FPR, TPR, p, iteration); 
  # auc_new <-  AUC_1(A, A_opt, FPR, TPR, p, iteration, auc);
  # return(list("RMSE_l0" = RMSE_Lasso_P, "ICOV_l0" = ICOV_P, "perf" =  perf, "AUC" = auc_new));
  return(list("RMSE_l0" = RMSE, "ICOV_l0" = ICOV, "RMSE_F" = RMSE_f, "ICOV_F" = ICOV_f));
}

RMSE_l0 <- function(true_cor, estimate_cor){
  p <- nrow(true_cor);
  Error <- abs(true_cor - estimate_cor);
  RMSE <- sum(sum(Error))/(p*(p-1));
  return(RMSE);
}

ICOV_l0 <- function(true_icov, estimate_icov){
  p <- nrow(true_icov);
  Error <- abs(true_icov - estimate_icov);
  RMSE <- sum(sum(Error))/(p*(p-1));
  return(RMSE);
}


RMSE_F <- function(true_cor, estimate_cor){
  p <- nrow(true_cor);
  Error <- abs(true_cor - estimate_cor)*abs(true_cor - estimate_cor);
  RMSE <- sqrt(sum(sum(Error)));
  return(RMSE);
}

ICOV_F <- function(true_icov, estimate_icov){
  p <- nrow(true_icov);
  Error <- abs(true_icov - estimate_icov)*abs(true_icov - estimate_icov);
  RMSE <- sqrt(sum(sum(Error)));
  return(RMSE);
}


ROCnew <- function(label, pred){
  
  pred <- prediction(pred, label)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred, 'auc')
  aucstd <- sd(simplify2array(auc@y.values))
  auc <- mean(simplify2array(auc@y.values))
  return(list("perf" = perf, "auc" = auc, "aucsd" = aucstd));
}


ROC <- function(true_adj, estimate_adj, FPR, TPR, p, iteration){
  
  A <- true_adj;
  A_opt <- estimate_adj;
  A_pred <- A_opt[upper.tri(A_opt,diag=FALSE)];
  A_label <- A[upper.tri(A,diag=FALSE)];
  pred <- prediction( A_pred, A_label);
  perf <- performance(pred, "tpr", "fpr");
  FPR[iteration+1] <- perf@x.values[[1]][2];
  TPR[iteration+1] <- perf@y.values[[1]][2];
  # TPR <- TPR[order(FPR, decreasing=FALSE)]
  # FPR <- sort(FPR, decreasing=F);
  perf@x.values[[1]] = FPR;
  perf@y.values[[1]] = TPR;
  return(perf);
}

AUC_1 <- function(true_adj, estimate_adj, FPR, TPR, p, iteration, auc_prev){
  A <- true_adj;
  A_opt <- estimate_adj;
  A_pred <- A_opt[upper.tri(A_opt,diag=FALSE)];
  A_label <- A[upper.tri(A,diag=FALSE)];
  pred <- prediction( A_pred, A_label);
  auc <- performance(pred, 'auc')@y.values[[1]];
  #auc_prev[iteration+1] <- performance(pred, 'auc')@y.values[[1]];
  #auc <- (auc_prev*(iteration-1) + auc)/iteration; 
  return(auc);
}

accEva <- function(label, pred){
  pred <- prediction(pred, label)
  perf <- performance(pred,"acc")
  acc <- matrix(0, length(label)) 
  for (i in 1:length(label)){
    index <- which.max(perf@y.values[[i]])
    acc[i] <- perf@y.values[[i]][index];
  }
  return(acc);
}

reproEva <- function(true_adj, test_adj){
  row_sum <- rowSums(true_adj);
  degree <- sort(row_sum, index.return = TRUE, decreasing=TRUE);
  degree_index <- degree$ix;
  num_taxa <- nrow(true_adj);
  percent = c(0.25, 0.5, 0.75, 1);
  reproError = matrix(0, 1, length(percent))
  for(i in 1:length(percent)){
    temp_index <- degree_index[1:round(num_taxa*percent[i])]
    #print(temp_index)
    num_edge_true <- sum(true_adj[temp_index,]) - round(num_taxa*percent[i]) 
    num_edge_test <- sum(test_adj[temp_index,]) - round(num_taxa*percent[i])
    error <- true_adj[temp_index,] -test_adj[temp_index,];
    #print(nrow(error))
    reproError[i] <- sum(abs(error))/round((num_taxa*num_taxa*percent[i]))
    #print(sum(abs(error)))
    #print(round(num_taxa*percent[i]))
  }

  return(reproError);
  
}
