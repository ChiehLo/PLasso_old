selection <- function(sample, type = "glasso_p", prior = 1, prior_type = TRUE, interaction_flag = FALSE){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample);
  
  lambda_min = 0.01;
  lambda_max = 1;
  interval = 0.01;
  numLambda = (lambda_max - lambda_min)/interval + 1;
  BIC_result = matrix(0, numLambda, 1);
  
  for (i in 1:numLambda){
    result <- glasso_p(sample, lambda_min = (lambda_min+(i-1)*interval), lambda_max = 10, prior, prior_type = prior_type, interaction_flag = interaction_flag);
    BIC_result[i] <- BIC(sample, result);
  }
  
  return(BIC_result);
  
  # A <- sample$adj;
  # A_opt <- result$adj;
  # 
  # perf <- ROC(A, A_opt, FPR, TPR, p, iteration); 
  # auc_new <-  AUC_1(A, A_opt, FPR, TPR, p, iteration, auc);
  # return(list("RMSE_l0" = RMSE_Lasso_P, "perf" =  perf, "AUC" = auc_new));
}

BIC <- function(sample, result){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample); 
  gamma <- 0.5;
  likelihood <- n/2*(log(det(result$icov)) - sum(diag(sample$var %*% result$icov)) );
  edge <- (sum(result$adj) - p)/2*log(n);
  extend <- 4*(sum(result$adj) - p)/2*log(p)*gamma;
  #print(likelihood)
  #print(edge)
  return(-2*likelihood + edge );
  #return(-2*likelihood);
}