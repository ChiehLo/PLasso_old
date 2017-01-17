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
library(ccrepe)
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
source("~/Desktop/HMP_project/my_R/CalculateCoOccurrence.R")
source('~/Desktop/HMP_project/my_R/real_data.R')
source("~/Desktop/HMP_project/REBACCA/REBACCA_main.R")
source("~/Desktop/HMP_project/REBACCA/REBACCA_simulation_model.R")
source('~/Desktop/HMP_project/CCLasso/CCLasso/R/cclasso.R')


run = 1;
eva = 0;
#load(file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/synthetic/cluster_100_200.RData");

glasso_p_indicator = TRUE;
cclasso_indicator = FALSE;
REBACCA_indicator = FALSE;
SparCC_indicator = FALSE;
mb_indicator = FALSE;
glasso_indicator = FALSE;
ccrepe_indicator = FALSE;

if (run ==1){
  OTUnum <- 100; # number of OTU
  Simulations <- 2e2; # number of sample
  Strength <- 0.2;
  numIt <- 100;
  sampleLabel <- matrix(0, numIt, OTUnum*OTUnum);
  #glasso_p
  prior <- 0.5;
  lambda_p <- 1; # glasso_p
  RMSEgp = matrix(0, numIt);
  ICOVgp = matrix(0, numIt);
  RMSEFgp = matrix(0, numIt);
  ICOVFgp = matrix(0, numIt);
  predgp = matrix(0, numIt, OTUnum*OTUnum);
  #labelgp = matrix(0, numIt, OTUnum*OTUnum);
  
  
  #cclasso
  
  RMSEcc = matrix(0, numIt);
  ICOVcc = matrix(0, numIt);
  RMSEFcc = matrix(0, numIt);
  ICOVFcc = matrix(0, numIt);
  predcc = matrix(0, numIt, OTUnum*OTUnum);
  #labelcc = matrix(0, numIt, OTUnum*OTUnum);
  
  
  #REBACCA
  
  RMSEr = matrix(0, numIt);
  ICOVr = matrix(0, numIt);
  RMSEFr = matrix(0, numIt);
  ICOVFr = matrix(0, numIt);
  predr = matrix(0, numIt, OTUnum*OTUnum);
  
  #SparCC
  
  RMSEs = matrix(0, numIt);
  ICOVs = matrix(0, numIt);
  RMSEFs = matrix(0, numIt);
  ICOVFs = matrix(0, numIt);
  preds = matrix(0, numIt, OTUnum*OTUnum);
  
  #SPIEC-EASI MB
  RMSEm = matrix(0, numIt);
  ICOVm = matrix(0, numIt);
  RMSEFm = matrix(0, numIt);
  ICOVFm = matrix(0, numIt);
  predm = matrix(0, numIt, OTUnum*OTUnum);
  
  
  #SPIEC-EASI glasso
  RMSEgl = matrix(0, numIt);
  ICOVgl = matrix(0, numIt);
  RMSEFgl = matrix(0, numIt);
  ICOVFgl = matrix(0, numIt);
  predgl = matrix(0, numIt, OTUnum*OTUnum);
  
  #ccrepe
  RMSEre = matrix(0, numIt);
  ICOVre = matrix(0, numIt);
  RMSEFre = matrix(0, numIt);
  ICOVFre = matrix(0, numIt);
  predre = matrix(0, numIt, OTUnum*OTUnum);
  
  ptm <- proc.time()
  for (i in 1:numIt){
    
    #sample generation
    sample <- graph_select(OTUnum = OTUnum, Simulations = Simulations, Strength = Strength, type = "AR4");
    sampleLabel[i,] <- abs(as.vector(sample$adj));
    print(i)
    if (glasso_p_indicator){
      select_lambda <- selection(sample, type="glasso_p", prior = prior);
      lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
      result_g <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "glasso_p");
      predgp[i,] <- abs(as.vector(result_g$icov));
      
      # calculate L1 distance
      evaluation_glasso_p <- evaluation(sample, result_g);
      RMSEgp[i] <- evaluation_glasso_p$RMSE_l0;
      ICOVgp[i] <- evaluation_glasso_p$ICOV_l0;
      RMSEFgp[i] <- evaluation_glasso_p$RMSE_F;
      ICOVFgp[i] <- evaluation_glasso_p$ICOV_F;
    }
    
    
    if (cclasso_indicator){
      result_c <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "cclasso");
      predcc[i,] <- abs(as.vector(result_c$icov));
      
      # calculate L1 distance
      evaluation_cclasso <- evaluation(sample, result_c);
      RMSEcc[i] <- evaluation_cclasso$RMSE_l0;
      ICOVcc[i] <- evaluation_cclasso$ICOV_l0;
      RMSEFcc[i] <- evaluation_cclasso$RMSE_F;
      ICOVFcc[i] <- evaluation_cclasso$ICOV_F;
    }
    
    if (REBACCA_indicator){
      result_r <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "REBACCA");
      predr[i,] <- abs(as.vector(result_r$adj));
      
      # calculate L1 distance
      evaluation_r <- evaluation(sample, result_r);
      RMSEr[i] <- evaluation_r$RMSE_l0;
      ICOVr[i] <- evaluation_r$ICOV_l0;
      RMSEFr[i] <- evaluation_r$RMSE_F;
      ICOVFr[i] <- evaluation_r$ICOV_F;
    }
    
    if (SparCC_indicator){
      result_s <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "sparcc");
      threshold <- 0.0;
      temp <- abs(result_s$cor) - threshold;
      temp[temp<0] <- 0;
      preds[i,] <- as.vector(temp);
      
      # calculate L1 distance
      evaluation_s <- evaluation(sample, result_s);
      RMSEs[i] <- evaluation_s$RMSE_l0;
      ICOVs[i] <- evaluation_s$ICOV_l0;
      RMSEFs[i] <- evaluation_s$RMSE_F;
      ICOVFs[i] <- evaluation_s$ICOV_F;
    }
    
    if(mb_indicator){
      result_m <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "SPIEC-EASI_mb");
      predm[i,] <- abs(as.vector(result_m$adj));
      
      # calculate L1 distance
      evaluation_m <- evaluation(sample, result_m);
      RMSEm[i] <- evaluation_m$RMSE_l0;
      ICOVm[i] <- evaluation_m$ICOV_l0;
      RMSEFm[i] <- evaluation_m$RMSE_F;
      ICOVFm[i] <- evaluation_m$ICOV_F;
    }
    
    if(glasso_indicator){
      result_gl <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "SPIEC-EASI_gl");
      predgl[i,] <- abs(as.vector(result_gl$icov));
      
      # calculate L1 distance
      evaluation_gl <- evaluation(sample, result_gl);
      RMSEgl[i] <- evaluation_gl$RMSE_l0;
      ICOVgl[i] <- evaluation_gl$ICOV_l0;
      RMSEFgl[i] <- evaluation_gl$RMSE_F;
      ICOVFgl[i] <- evaluation_gl$ICOV_F;
    }
    
    if(ccrepe_indicator){
      result_re <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "ccrepe");
      predre[i,] <-  abs(as.vector(result_re$cov));
      # calculate L1 distance
      evaluation_re <- evaluation(sample, result_re);
      RMSEre[i] <- evaluation_re$RMSE_l0;
      ICOVre[i] <- evaluation_re$ICOV_l0;
      RMSEFre[i] <- evaluation_re$RMSE_F;
      ICOVFre[i] <- evaluation_re$ICOV_F;
    }
    
  }
  print(proc.time() - ptm)
}

if(eva==1){
  if (glasso_p_indicator){
    templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
    temppred <- split(t(predgp), rep(1:nrow(predgp), each = ncol(predgp)))
    glasso_p_ROC <- ROCnew(templabel, temppred)
    perfgp <- glasso_p_ROC$perf;
    aucgp <- glasso_p_ROC$auc;
    aucstdgp <- glasso_p_ROC$aucsd
    print("glasso_p:")
    #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(ICOVgp), sd(ICOVgp), mean(ICOVFgp), sd(ICOVFgp), aucgp, aucstdgp))
    print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEgp), sd(RMSEgp), mean(RMSEFgp), sd(RMSEFgp), aucgp, aucstdgp))
    #plot(perfgp,col="grey82",lty=3)
    plot(perfgp,lwd=3,avg="threshold", col="red", yaxis.at=c(-1, 2.0), xaxis.at = c(-1, 2),ann = FALSE)
    #plot(perfgp,lwd=3,avg="threshold", col="red", yaxis.at=c(-1, 2.0), xaxis.at = c(-1, 2),ann = FALSE)
    axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), NA, cex.axis=.7, font=1, tck=.03)
    axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), NA, cex.axis=.7, font=1, tck=.03)
    #save(perfgp, file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/temp/test.RData")
    #load(file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/temp/test.RData")
  }
  if (cclasso_indicator){
    templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
    temppred <- split(t(predcc), rep(1:nrow(predcc), each = ncol(predcc)))
    cclasso_ROC <- ROCnew(templabel, temppred)
    perfcc <- cclasso_ROC$perf;
    auccc <- cclasso_ROC$auc;
    aucstdcc <- cclasso_ROC$aucsd
    print("cclasso:")
    #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEcc), sd(RMSEcc), mean(RMSEFcc), sd(RMSEFcc), auccc, aucstdcc))
    #plot(perfcc,col="grey82",lty=3)
    plot(perfcc,lwd=3,avg="threshold",add=TRUE, col = "darkgreen", lty=2)
  }
  
  if (REBACCA_indicator){
    templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
    temppred <- split(t(predr), rep(1:nrow(predr), each = ncol(predr)))
    REBACCA_ROC <- ROCnew(templabel, temppred)
    perfr <- REBACCA_ROC$perf;
    aucr <- REBACCA_ROC$auc;
    aucstdr <- REBACCA_ROC$aucsd
    print("REBACCA:")
    #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEr), sd(RMSEr), mean(RMSEFr), sd(RMSEFr), aucr, aucstdr))
    #plot(perfr,col="grey82",lty=3)
    plot(perfr,lwd=3,avg="threshold",add=TRUE, col = "gold2", lty=2)
  }
  
  if (SparCC_indicator){
    templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
    temppred <- split(t(preds), rep(1:nrow(preds), each = ncol(preds)))
    SparCC_ROC <- ROCnew(templabel, temppred)
    perfs <- SparCC_ROC$perf;
    aucs <- SparCC_ROC$auc;
    aucstds <- SparCC_ROC$aucsd
    print("SparCC:")
    #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEs), sd(RMSEs), mean(RMSEFs), sd(RMSEFs), aucs, aucstds))
    #plot(perfs,col="grey82",lty=3)
    plot(perfs,lwd=3,avg="threshold",add=TRUE, col = "blue", lty=2)
  }
  
  if (mb_indicator){
    templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
    temppred <- split(t(predm), rep(1:nrow(predm), each = ncol(predm)))
    mb_ROC <- ROCnew(templabel, temppred)
    perfmb <- mb_ROC$perf;
    aucmb <- mb_ROC$auc;
    aucstdmb <- mb_ROC$aucsd
    print("mb:")
    #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(ICOVm), sd(ICOVm), mean(ICOVFm), sd(ICOVFm), aucmb, aucstdmb))
    #plot(perfmb,col="grey82",lty=3)
    plot(perfmb,lwd=3,avg="threshold",add=TRUE, col = "cyan3", lty=2)
  }
  
  if (glasso_indicator){
    templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
    temppred <- split(t(predgl), rep(1:nrow(predgl), each = ncol(predgl)))
    gl_ROC <- ROCnew(templabel, temppred)
    perfgl <- gl_ROC$perf;
    aucgl <- gl_ROC$auc;
    aucstdgl <- gl_ROC$aucsd
    print("glasso:")
    #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(ICOVgl), sd(ICOVgl), mean(ICOVFgl), sd(ICOVFgl), aucgl, aucstdgl))
    #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEgl), sd(RMSEgl), mean(RMSEFgl), sd(RMSEFgl), aucgl, aucstdgl))
    #plot(perfgl,col="grey82",lty=3)
    plot(perfgl,lwd=3,avg="threshold",add=TRUE, col = "gray8", lty=2)
  }
  
  if (ccrepe_indicator){
    templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
    temppred <- split(t(predre), rep(1:nrow(predre), each = ncol(predre)))
    re_ROC <- ROCnew(templabel, temppred)
    perfre <- re_ROC$perf;
    aucre <- re_ROC$auc;
    aucstdre <- re_ROC$aucsd
    print("ccrepe:")
    #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(ICOVre), sd(ICOVre), mean(ICOVFre), sd(ICOVFre), aucre, aucstdre))
    #plot(perfre,col="grey82",lty=3)
    plot(perfre,lwd=3,avg="threshold",add=TRUE, col = "purple", lty=2)
  }
  #save(sampleLabel, predcc, predgl, predgp, predm, predr, predre, preds, perfcc, perfgl, perfgp, perfmb, perfr, perfre, perfs, file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/synthetic/scale_free_100_200.RData")
}

#plot(perf,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE)

# TPRgp <- TPRgp[order(FPRgp, decreasing=FALSE)]
# FPRgp <- sort(FPRgp, decreasing=F);
# evaluation_glasso_p$perf@x.values[[1]] = FPRgp;
# evaluation_glasso_p$perf@y.values[[1]] = TPRgp;
#plot(evaluation_glasso_p$perf)


