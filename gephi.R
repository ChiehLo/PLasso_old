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
HMMCP = FALSE;
HMQCP = FALSE;
Synthetic = TRUE;

if (Synthetic == TRUE){
  sample <- graph_select(OTUnum = 20, Simulations = 20, Strength = Strength, type = "scale-free");
  test<- graph_from_adjacency_matrix(sample$adj, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
  gexf2 <- igraph.to.gexf(test)
  #print(gexf2, "/Users/chiehlo/Desktop/HMP_project/my_R/graph/true_10.gexf", replace=T)
  Adj <- sample$adj
  for(i in 1:10){
    for (j in i:10){
      if(Adj[i,j]==1 && runif(1, min=0, max=1) > 0.5 && i!=j){
        Adj[i,j] <- 0;
        Adj[j,i] <- 0;
      }
    }
  }
  test<- graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
  gexf2 <- igraph.to.gexf(test)
  #print(gexf2, "/Users/chiehlo/Desktop/HMP_project/my_R/graph/prior_10.gexf", replace=T)
}


## HMP
if (HMASM == TRUE){
  pre_path = "/Users/chiehlo/Desktop/HMP_project/data/HMASM/";
  bodysite = "Tongue_dorsum"
  label_file = "/output/SpeciesName1.txt";
  label_path = paste(pre_path, bodysite, label_file, sep= ''); 
  Adj <- result_g$adj;
  # label 
  #label_temp <- t(as.matrix(read.table("/Users/chiehlo/Desktop/HMP_project/data/HMASM/Anterior_nares/output/SpeciesName1.txt")));
  label_temp <- t(as.matrix(read.table(label_path)));
  # label <- c();
  # for (i in 1:ncol(label_temp)){
  #   label[i] <-  paste(label_temp[1:2, i], collapse = "_")
  # }
  label <- label_temp
  rownames(Adj) <- label;
  colnames(Adj) <- label;
  test<- graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
  gexf2 <- igraph.to.gexf(test)
  output_path = paste("/Users/chiehlo/Desktop/HMP_project/my_R/graph/", bodySite, "_HMASM.gexf", sep= ''); 
  print(gexf2, output_path, replace=T)
  #print(gexf2, "/Users/chiehlo/Desktop/HMP_project/my_R/graph/Anterior_nares_HMASM.gexf", replace=T)
}

if (HMMCP == TRUE){
  Adj <- result_g$adj;
  # label
  pre_path = "/Users/chiehlo/Desktop/HMP_project/data/HMMCP/v35/";
  bodysite = "Tongue_dorsum"
  label_file = "/genus_name_filter.txt";
  label_path = paste(pre_path, bodysite, label_file, sep= ''); 
  label_temp <- t(as.matrix(read.table(label_path)));
  # label_temp <- t(as.matrix(read.table("/Users/chiehlo/Desktop/HMP_project/data/HMMCP/v35/Anterior_nares/genus_name_filter.txt")));
  label <- label_temp;

  rownames(Adj) <- label;
  colnames(Adj) <- label;
  test<- graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
  gexf2 <- igraph.to.gexf(test)
  output_path = paste("/Users/chiehlo/Desktop/HMP_project/my_R/graph/", bodySite, "_HMMCP.gexf", sep= ''); 
  print(gexf2, output_path, replace=T)
  # print(gexf2, "/Users/chiehlo/Desktop/HMP_project/my_R/graph/Anterior_nares_HMMCP.gexf", replace=T)
}

if (HMQCP == TRUE){
  Adj <- result_g$adj;
  # label 
  pre_path = "/Users/chiehlo/Desktop/HMP_project/data/HMQCP/v35/";
  bodysite = "Tongue_dorsum"
  label_file = "/genus_name_filter.txt";
  label_path = paste(pre_path, bodysite, label_file, sep= ''); 
  label_temp <- t(as.matrix(read.table(label_path)));
  # label_temp <- t(as.matrix(read.table("/Users/chiehlo/Desktop/HMP_project/data/HMQCP/v35/Anterior_nares/genus_name_filter.txt")));
  label <- label_temp;
  
  rownames(Adj) <- label;
  colnames(Adj) <- label;
  test<- graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
  gexf2 <- igraph.to.gexf(test)
  output_path = paste("/Users/chiehlo/Desktop/HMP_project/my_R/graph/", bodySite, "_HMQCP.gexf", sep= ''); 
  print(gexf2, output_path, replace=T)
  # print(gexf2, "/Users/chiehlo/Desktop/HMP_project/my_R/graph/Anterior_nares_HMQCP.gexf", replace=T)
}
