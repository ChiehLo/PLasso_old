library(MASS)
library(Matrix)
library(matrixcalc)
library(glasso)
library(ROCR)
library(huge)
library(MCMCpack)
library(glmnet)
library(parallel)
library(R.matlab)
source('~/Desktop/HMP_project/SpiecEasi-master/R/community_graph.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/fitdistr.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/mvdistributions.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/normalization.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/plotNet.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/roc.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/spaRcc.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/SparseICov.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/utilities.R')
source('~/Desktop/HMP_project/my_R/syn_data_generate.R');
source('~/Desktop/HMP_project/my_R/algorithm.R');
source('~/Desktop/HMP_project/my_R/evaluation.R');
source('~/Desktop/HMP_project/my_R/selection.R');
source('~/Desktop/HMP_project/CCLasso/CCLasso/R/cclasso.R')
source('~/Desktop/HMP_project/my_R/real_data.R')
source("/Users/chiehlo/Desktop/HMP_project/REBACCA/REBACCA_main.R")
source("/Users/chiehlo/Desktop/HMP_project/REBACCA/REBACCA_simulation_model.R")
source("/Users/chiehlo/Desktop/HMP_project/my_R/CalculateCoOccurrence.R")

OTUnum <- 50; # number of OTU
Simulations <- 5e2; # number of sample
Strength <- 0.2;
numIt <- 1;
FPR = matrix(0, numIt+2,1);
TPR = matrix(0, numIt+2,1);
FPR[1] = 0;
TPR[1] = 0;
FPR[numIt+2] = 1;
TPR[numIt+2] = 1;
AUC <- 0;
prior <- 0.1;
lambda_p <- 1; # glasso_p

for (i in 1:numIt){


  sample <- graph_select(OTUnum = OTUnum, Simulations = Simulations, Strength = Strength, type = "random");
  select_lambda <- selection(sample, type="glasso_p", prior = prior);
  lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
  result_c <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "cclasso");
  result_g <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "glasso_p");
  #result_r <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "REBACCA");


  #result_s <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "sparcc");



  #evaluation
  evaluation_glasso_p <- evaluation(sample, result_g, i, FPR, TPR, AUC);
  FPR <- evaluation_glasso_p$perf@x.values[[1]];
  TPR <- evaluation_glasso_p$perf@y.values[[1]];
  AUC <- evaluation_glasso_p$AUC;

  evaluation_cclasso <- evaluation(sample, result_c, i, FPR, TPR, AUC);
  #evaluation_r <- evaluation(sample, result_r, i, FPR, TPR, AUC);
}
TPR <- TPR[order(FPR, decreasing=FALSE)]
FPR <- sort(FPR, decreasing=F);
evaluation_glasso_p$perf@x.values[[1]] = FPR;
evaluation_glasso_p$perf@y.values[[1]] = TPR;
plot(evaluation_glasso_p$perf)

# for (i in 1:numIt){
#   #association <- occurance();
#   interaction <- interaction();
#   #sample <- readOTU(level = "genus", site = "stool", 1 - No_association);
#   sample <- readOTU(level = "species", site = "stool", interaction);
#   select_lambda <- selection(sample, type="glasso_p", prior = prior, prior_type = FALSE);
#   lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
#   #result_c <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "cclasso");
#   result_g <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "glasso_p", prior_type = FALSE);
#   #result_r <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "REBACCA");
#   
# 
#   #result_s <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "sparcc");
# 
# 
# 
#   #evaluation
#   # evaluation_glasso_p <- evaluation(sample, result_g, i, FPR, TPR, AUC);
#   # FPR <- evaluation_glasso_p$perf@x.values[[1]];
#   # TPR <- evaluation_glasso_p$perf@y.values[[1]];
#   # AUC <- evaluation_glasso_p$AUC;
# 
#   #evaluation_cclasso <- evaluation(sample, result_c, i, FPR, TPR, AUC);
#   #evaluation_r <- evaluation(sample, result_r, i, FPR, TPR, AUC);
# }
# # TPR <- TPR[order(FPR, decreasing=FALSE)]
# # FPR <- sort(FPR, decreasing=F);
# # evaluation_glasso_p$perf@x.values[[1]] = FPR;
# # evaluation_glasso_p$perf@y.values[[1]] = TPR;
# # plot(evaluation_glasso_p$perf)
se.est <- spiec.easi(sample$sample, method='mb', lambda.min.ratio=1e-2, nlambda=15)
elist.mb <- symBeta(getOptBeta(se.est), mode='maxabs')
result <- as.matrix(elist.mb)
adj <-as.matrix(se.est$refit)