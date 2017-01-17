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
library(matrixStats)
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
source('~/Desktop/HMP_project/CCLasso/CCLasso/R/SparCC.R')
source('~/Desktop/HMP_project/my_R/real_data.R')
source("/Users/chiehlo/Desktop/HMP_project/REBACCA/REBACCA_main.R")
source("/Users/chiehlo/Desktop/HMP_project/REBACCA/REBACCA_simulation_model.R")
source("/Users/chiehlo/Desktop/HMP_project/my_R/CalculateCoOccurrence.R")

bodysite <- c('Anterior_nares', 'Buccal_mucosa', 'Stool', 'Supragingival_plaque', 'Tongue_dorsum')
bodysite <- c('Tongue_dorsum')
HMASM = TRUE;
HMMCP = TRUE;
HMQCP = TRUE;
repro = FALSE;
if (HMASM == TRUE){
  for (i in 1:length(bodysite)){
    bodySite = bodysite[i];
    pre_path = "/Users/chiehlo/Desktop/HMP_project/data/HMASM/";
    interaction_file = "/output/species_interaction.csv";
    interaction_path = paste(pre_path, bodySite, interaction_file, sep= ''); 
    #interaction_path =  "/Users/chiehlo/Desktop/HMP_project/data/HMASM/Buccal_mucosa/output/species_interaction.csv"
    interaction_matrix <- interaction(interaction_path);
    occu_path = paste(pre_path, bodySite, "/output/association.csv", sep = '');
    #occu_path = "/Users/chiehlo/Desktop/HMP_project/data/HMASM/Buccal_mucosa/output/association.csv"
    No_association <- occurance(occu_path);
    # #sample <- readOTU(level = "genus", site = "stool", 1 - No_association);
    count_path = paste(pre_path, bodySite, "/output/SpeciesCount.txt", sep = '');
    #count_path = "/Users/chiehlo/Desktop/HMP_project/data/HMASM/Buccal_mucosa/output/SpeciesCount.txt"
    sample <- readOTU(level = "species", count_path = count_path, 1-No_association, interaction_matrix);
    prior = 0;
    interaction_flag = TRUE;
    print(interaction_flag)
    select_lambda <- selection(sample, type="glasso_p", prior = prior, interaction_flag = interaction_flag);
    lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
    result_g <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "glasso_p", interaction_flag = interaction_flag);
    n_boot = 20;
    result_c <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "cclasso", n_boot = n_boot, counts = F);
    nInt <- 1;
    reproduceErrorArray_g = matrix(0,nInt, 4);
    reproduceErrorArray_c = matrix(0,nInt, 4);
    if (repro == TRUE){
      for (k in 1:nInt){
        print(k)
        num_sample <- nrow(sample$sample)
        subsample <- sample;
        sample_index <- sort(sample.int(num_sample, round(0.5*num_sample)))
        subsample$sample <- subsample$sample[sample_index,]
        lasso_sub <- fraction_data_process(subsample$sample, 1 - No_association, interaction_matrix)
        select_lambda <- selection(lasso_sub, type="glasso_p", prior = prior, interaction_flag = interaction_flag);
        lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
        result_g_sample <- algorithm_select(lasso_sub, lambda_min = lambda_p, prior = prior, type = "glasso_p", interaction_flag = interaction_flag);
        result_c_sample <- algorithm_select(subsample, lambda_min = lambda_p, prior = prior, type = "cclasso", n_boot = n_boot, counts = F);
        #result_s_sample <- algorithm_select(subsample, lambda_min = lambda_p, prior = prior, type = "sparcc");
        reproduceError_g <- reproEva(result_g$adj, result_g_sample$adj)
        reproduceError_c <- reproEva(result_c$adj, result_c_sample$adj)
        #reproduceError_s <- reproEva(result_s$adj, result_s_sample$adj)
        reproduceErrorArray_g[k,] <- reproduceError_g;
        reproduceErrorArray_c[k,] <- reproduceError_c;
        #print(reproduceError_g)
        #print(reproduceError_c)
        #print(reproduceError_s)        
      }
      Mean_g <- colMeans(1-reproduceErrorArray_g)
      Std_g <- colSds(1-reproduceErrorArray_g)
      Mean_c <- colMeans(1-reproduceErrorArray_c)
      Std_c <- colSds(1-reproduceErrorArray_c)
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", Mean_g[1], Std_g[1], Mean_c[1], Std_c[1], Mean_g[2], Std_g[2] , Mean_c[2], Std_c[2], Mean_g[3], Std_g[3],Mean_c[3], Std_c[3], Mean_g[4], Std_g[4],Mean_c[4], Std_c[4]))
      #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", Mean_c[1], Std_c[1], Mean_c[2], Std_c[2] , Mean_c[3], Std_c[3], Mean_c[4], Std_c[4]))
    }
    
    # result_c <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "cclasso");
    # result_s <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "sparcc");
  }  
}

if (HMMCP == TRUE){
  for (i in 1:length(bodysite)){
    bodySite = bodysite[i];
    pre_path = "/Users/chiehlo/Desktop/HMP_project/data/HMMCP/v35/";
    occu_path = paste(pre_path, bodySite, "/association.csv", sep = '');
    No_association <- occurance(occu_path);
    
    
    count_path = paste(pre_path, bodySite, "/genus_count_filter.csv", sep = '');
    sample <- readOTU(level = "genus", count_path = count_path, 1 - No_association);
    #sample <- readOTU(level = "species", count_path = count_path, 1-No_association, interaction_matrix);
    prior = 0;
    interaction_flag = FALSE;
    print(interaction_flag)
    select_lambda <- selection(sample, type="glasso_p", prior = prior, interaction_flag = interaction_flag);
    lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
    result_g <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "glasso_p", interaction_flag = interaction_flag);
    result_g_MM <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "glasso_p", interaction_flag = interaction_flag);
    #result_c <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "cclasso");
    #result_s <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "sparcc");
    nInt <- 20;
    reproduceErrorArray_g = matrix(0,nInt, 4);
    reproduceErrorArray_c = matrix(0,nInt, 4);
    n_boot = 20;
    if (repro == TRUE){
      for (k in 1:nInt){
        print(k)
        num_sample <- nrow(sample$sample)
        subsample <- sample;
        sample_index <- sort(sample.int(num_sample, round(0.5*num_sample)))
        subsample$sample <- subsample$sample[sample_index,]
        lasso_sub <- count_data_process(subsample$sample, 1 - No_association)
        select_lambda <- selection(lasso_sub, type="glasso_p", prior = prior, interaction_flag = interaction_flag);
        lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
        result_g_sample <- algorithm_select(lasso_sub, lambda_min = lambda_p, prior = prior, type = "glasso_p", interaction_flag = interaction_flag);
        result_c_sample <- algorithm_select(subsample, lambda_min = lambda_p, prior = prior, type = "cclasso", n_boot = n_boot, counts = T);
        #result_s_sample <- algorithm_select(subsample, lambda_min = lambda_p, prior = prior, type = "sparcc");
        reproduceError_g <- reproEva(result_g$adj, result_g_sample$adj)
        reproduceError_c <- reproEva(result_c$adj, result_c_sample$adj)
        #reproduceError_s <- reproEva(result_s$adj, result_s_sample$adj)
        reproduceErrorArray_g[k,] <- reproduceError_g;
        reproduceErrorArray_c[k,] <- reproduceError_c;
        #print(reproduceError_g)
        #print(reproduceError_c)
        #print(reproduceError_s)        
      }
      Mean_g <- colMeans(1-reproduceErrorArray_g)
      Std_g <- colSds(1-reproduceErrorArray_g)
      Mean_c <- colMeans(1-reproduceErrorArray_c)
      Std_c <- colSds(1-reproduceErrorArray_c)
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", Mean_g[1], Std_g[1], Mean_c[1], Std_c[1], Mean_g[2], Std_g[2] , Mean_c[2], Std_c[2], Mean_g[3], Std_g[3],Mean_c[3], Std_c[3], Mean_g[4], Std_g[4],Mean_c[4], Std_c[4]))
      #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", Mean_c[1], Std_c[1], Mean_c[2], Std_c[2] , Mean_c[3], Std_c[3], Mean_c[4], Std_c[4]))
    }
    
  }  
}

if (HMQCP == TRUE){
  for (i in 1:length(bodysite)){
    bodySite = bodysite[i];
    pre_path = "/Users/chiehlo/Desktop/HMP_project/data/HMQCP/v35/";
    occu_path = paste(pre_path, bodySite, "/association.csv", sep = '');
    No_association <- occurance(occu_path);
    
    
    count_path = paste(pre_path, bodySite, "/genus_count_filter.csv", sep = '');
    sample <- readOTU(level = "genus", count_path = count_path, 1 - No_association);
    #sample <- readOTU(level = "species", count_path = count_path, 1-No_association, interaction_matrix);
    prior = 0;
    interaction_flag = FALSE;
    print(interaction_flag)
    select_lambda <- selection(sample, type="glasso_p", prior = prior, interaction_flag = interaction_flag);
    lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
    result_g <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "glasso_p", interaction_flag = interaction_flag);
    result_g_MQ <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "glasso_p", interaction_flag = interaction_flag);
    #result_c <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "cclasso", counts = T);
    # result_s <- algorithm_select(sample, lambda_min = lambda_p, prior = prior, type = "sparcc");
    nInt <- 20;
    reproduceErrorArray_g = matrix(0,nInt, 4);
    reproduceErrorArray_c = matrix(0,nInt, 4);
    n_boot = 20;
    if (repro == TRUE){
      for (k in 1:nInt){
        print(k)
        num_sample <- nrow(sample$sample)
        subsample <- sample;
        sample_index <- sort(sample.int(num_sample, round(0.5*num_sample)))
        subsample$sample <- subsample$sample[sample_index,]
        lasso_sub <- count_data_process(subsample$sample, 1 - No_association)
        select_lambda <- selection(lasso_sub, type="glasso_p", prior = prior, interaction_flag = interaction_flag);
        lambda_p <- (which.min(select_lambda)-1)*0.01 + 0.01;
        result_g_sample <- algorithm_select(lasso_sub, lambda_min = lambda_p, prior = prior, type = "glasso_p", interaction_flag = interaction_flag);
        result_c_sample <- algorithm_select(subsample, lambda_min = lambda_p, prior = prior, type = "cclasso", n_boot = n_boot, counts = F);
        #result_s_sample <- algorithm_select(subsample, lambda_min = lambda_p, prior = prior, type = "sparcc");
        reproduceError_g <- reproEva(result_g$adj, result_g_sample$adj)
        reproduceError_c <- reproEva(result_c$adj, result_c_sample$adj)
        #reproduceError_s <- reproEva(result_s$adj, result_s_sample$adj)
        reproduceErrorArray_g[k,] <- reproduceError_g;
        reproduceErrorArray_c[k,] <- reproduceError_c;
        #print(reproduceError_g)
        #print(reproduceError_c)
        #print(reproduceError_s)        
      }
      Mean_g <- colMeans(1-reproduceErrorArray_g)
      Std_g <- colSds(1-reproduceErrorArray_g)
      Mean_c <- colMeans(1-reproduceErrorArray_c)
      Std_c <- colSds(1-reproduceErrorArray_c)
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", Mean_g[1], Std_g[1], Mean_c[1], Std_c[1], Mean_g[2], Std_g[2] , Mean_c[2], Std_c[2], Mean_g[3], Std_g[3],Mean_c[3], Std_c[3], Mean_g[4], Std_g[4],Mean_c[4], Std_c[4]))
      #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", Mean_c[1], Std_c[1], Mean_c[2], Std_c[2] , Mean_c[3], Std_c[3], Mean_c[4], Std_c[4]))
    }
    
  }  
}
