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
library(igraph)
library(rgexf)
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

HMASM = FALSE;
HMMCP = TRUE;
HMQCP = TRUE;
Synthetic = TRUE;
evaluate_ranking = TRUE;


## HMP


if (HMMCP == TRUE){
  Adj_MM <- result_g_MM$adj;
  # label
  pre_path = "/Users/chiehlo/Desktop/HMP_project/data/HMMCP/v35/";
  bodysite = "Tongue_dorsum"
  label_file = "/genus_name_filter.txt";
  label_path = paste(pre_path, bodysite, label_file, sep= ''); 
  label_temp <- t(as.matrix(read.table(label_path)));
  # label_temp <- t(as.matrix(read.table("/Users/chiehlo/Desktop/HMP_project/data/HMMCP/v35/Anterior_nares/genus_name_filter.txt")));
  label <- label_temp;
  
  #rownames(Adj) <- label;
  #colnames(Adj) <- label;
  label_MM <- label_temp;
}

if (HMQCP == TRUE){
  Adj_MQ <- result_g_MQ$adj;
  # label 
  pre_path = "/Users/chiehlo/Desktop/HMP_project/data/HMQCP/v35/";
  bodysite = "Tongue_dorsum"
  label_file = "/genus_name_filter.txt";
  label_path = paste(pre_path, bodysite, label_file, sep= ''); 
  label_temp <- t(as.matrix(read.table(label_path)));
  # label_temp <- t(as.matrix(read.table("/Users/chiehlo/Desktop/HMP_project/data/HMQCP/v35/Anterior_nares/genus_name_filter.txt")));
  label <- label_temp;
  
  #rownames(Adj) <- label;
  #colnames(Adj) <- label;
  label_MQ <- label_temp;
}

if(evaluate_ranking){
  label_all<- c(label_MM, label_MQ);
  label_all <- unique(label_all);
  ranking_matrix = matrix(0, length(label_all), 2);
  for(i in 1:length(label_all)){
    MM_index = which(label_all[i] == label_MM);
    MQ_index = which(label_all[i] == label_MQ);
    if(length(MM_index)>0){
      ranking_matrix[i,1] = sum(Adj_MM[MM_index,])
    }
    else{
      ranking_matrix[i,1] = NA; 
    }
    if(length(MQ_index)>0){
      ranking_matrix[i,2] = sum(Adj_MQ[MQ_index,])
    }
    else{
      ranking_matrix[i,2] = NA;
    }
  }
  ranking_matrix[,1]<-rank(-ranking_matrix[,1], na.last = "keep")
  ranking_matrix[,2]<-rank(-ranking_matrix[,2], na.last = "keep")
  print(cor(ranking_matrix, use="complete.obs", method="spearman")) 
  #print(ranking_matrix)
}




