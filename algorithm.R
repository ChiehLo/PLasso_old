algorithm_select <- function(sample, lambda_min = 0.01, lambda_max = 10, prior = 1, type = "glasso_p", prior_type = TRUE, interaction_flag = FALSE, n_boot = 20, counts = T){
  if(type=="glasso_p"){
    result <- glasso_p(sample, lambda_min, lambda_max, prior, prior_type = prior_type, interaction_flag = interaction_flag)
  }
  else if(type == "sparcc"){
    result <- sparcc_algo(sample);
  }
  else if(type == "cclasso"){
    result <- cclasso_algo(sample, counts = counts, n_boot = n_boot);
  }
  else if(type == "REBACCA"){
    result <- REBACCA_algo(sample);
  }
  else if(type == "SPIEC-EASI_mb"){
    result <- SPIEC_mb(sample);
  }
  else if(type == "SPIEC-EASI_gl"){
    result <- SPIEC_gl(sample);
  }
  else if(type == "ccrepe"){
    result <- ccrepe_algo(sample)
  }
  return(result);
}


glasso_p <- function(sample, lambda_min = 0.01, lambda_max = 10, prior = 1, prior_type, interaction_flag){
  s <- sample$var;
  adj <- sample$adj;
  if(interaction_flag==TRUE){
    interaction <- sample$inter;
  }
  p <- nrow(s);
  print(lambda_min);
  rho <- lambda_min*matrix(1, p, p); #gLasso regularizer
  if (prior_type == TRUE){
    rho <- lambda_min*matrix(1, p, p); #gLasso regularizer
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 0 && runif(1, min=0, max=1) > prior){
          rho[i,j] <- lambda_max;
          rho[j,i] <- lambda_max;
        }
      }
    }
  }
  if(interaction_flag == TRUE){
    for (i in 1:p){
      for (j in 1:p){
        if (abs(interaction[i,j]) == 1 && runif(1, min=0, max=1) > prior && i!=j){
          rho[i,j] <- 0.00001;
          rho[j,i] <- 0.00001;
        }
      }
    }
  }
  #print(rho)
  a<-glasso(s, rho=rho)
  C_opt <- a$w;
  CI_opt <- a$wi;
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(C_opt[i,i]);
  }
  D <- diag(d);
  R_opt <- D%*%C_opt%*%D;
  
  A_opt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(CI_opt[i,j]) > 0){
        A_opt[i,j] <- 1;
        A_opt[j,i] <- 1;
      }
    }
  }
  
  
  return(list("cor" = R_opt, "cov" = C_opt, "icov" = CI_opt, "adj" = A_opt ));
}

sparcc_algo <- function(sample){
  s <- sample$var;
  adj <- sample$adj;
  p <- nrow(s);
  
  sparcc_result <- sparcc(sample$sample) #need to put count data
  R_opt <- sparcc_result$Cor;
  C_opt <- sparcc_result$Cov;
  CI_opt <- ginv(C_opt);
  #sparcc_R <- sparcc_result$Cor
  A_opt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(R_opt[i,j]) > 0.01){
        A_opt[i,j] <- 1;
        A_opt[j,i] <- 1;
      }
    }
  }
  return(list("cor" = R_opt, "cov" = C_opt, "icov" = CI_opt, "adj" = A_opt ));
  #return(sparcc_result);
}

cclasso_algo <- function(sample, counts, n_boot){
  s <- sample$var;
  adj <- sample$adj;
  p <- nrow(s);
  if(counts){
    cclasso_result <- cclasso(sample$sample, counts = counts, n_boot = n_boot) #need to put count data 
  }
  else{
    cclasso_result <- cclasso(sample$fraction, counts = counts, n_boot = n_boot) #need to put count data
  }
  R_opt <- cclasso_result$cor_w;
  d <- 0;
  for(i in 1:p){
    d[i] <- sqrt(cclasso_result$var_w[i]);
  }
  D <- diag(d);
  C_opt <- D%*%R_opt%*%D ;
  CI_opt <- ginv(C_opt);
  A_opt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(R_opt[i,j]) > 0.1){
        A_opt[i,j] <- 1;
        A_opt[j,i] <- 1;
      }
    }
  }
  return(list("cor" = R_opt, "cov" = C_opt, "icov" = CI_opt, "adj" = A_opt ));
  #return(sparcc_result);
}


REBACCA_algo <- function(sample){
  s <- sample$var;
  adj <- sample$adj;
  p <- nrow(s);
  
  x.rslt = rebacca(t(sample$frac), nbootstrap=50, N.cores=4)
  
  # estimate correlation
  tau = stability_cutoff(x.rslt$Stability, x.rslt$q, B=50, FWER=0.05)
  x.adj = sscore2adjmatrix(x.rslt$Stability, tau)
  diag(x.adj) <- 1;
  x.est = rebacca_adjm2corr(t(sample$frac), x.adj)
  
  R_opt <- x.est$corr;
  C_opt <- x.est$cov;
  CI_opt <- ginv(C_opt);
  A_opt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(x.adj[i,j]) > 0){
        A_opt[i,j] <- 1;
        A_opt[j,i] <- 1;
      }
    }
  }
  return(list("cor" = R_opt, "cov" = C_opt, "icov" = CI_opt, "adj" = A_opt ));
  #return(sparcc_result);
}

SPIEC_mb <- function(sample, lambda_min = 0.01, lambda_max = 10, prior = 1, prior_type){
  s <- sample$var;
  adj <- sample$adj;
  p <- nrow(s);

  se.est <- spiec.easi(sample$sample, method='mb', lambda.min.ratio=1e-2, nlambda=15)
  elist.mb <- symBeta(getOptBeta(se.est), mode='maxabs')
  result <- as.matrix(elist.mb)
  adj <-as.matrix(se.est$refit)

  C_opt <- matrix(0,p, p);
  CI_opt <- matrix(0, p, p);

  R_opt <- matrix(0, p, p);
  R_opt <- result;
  A_opt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(adj[i,j]) > 0){
        A_opt[i,j] <- 1;
        A_opt[j,i] <- 1;
      }
    }
  }
  diag(A_opt) <- 1;
  
  return(list("cor" = R_opt, "cov" = C_opt, "icov" = CI_opt, "adj" = A_opt ));
}

SPIEC_gl <- function(sample, lambda_min = 0.01, lambda_max = 10, prior = 1, prior_type){
  s <- sample$var;
  adj <- sample$adj;
  p <- nrow(s);
  
  se.est <- spiec.easi(sample$sample, method='glasso', lambda.min.ratio=1e-2, nlambda=15)
  print(se.est$opt.lambda)
  C_opt <- as.matrix(se.est$opt.cov);
  CI_opt <- as.matrix(se.est$opt.icov);
  
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(C_opt[i,i]);
  }
  D <- diag(d);
  R_opt <- D%*%C_opt%*%D;
  
  A_opt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(CI_opt[i,j]) > 0){
        A_opt[i,j] <- 1;
        A_opt[j,i] <- 1;
      }
    }
  }
  diag(A_opt) <- 1;
  
  return(list("cor" = R_opt, "cov" = C_opt, "icov" = CI_opt, "adj" = A_opt ));
}

ccrepe_algo <- function(sample, lambda_min = 0.01, lambda_max = 10, prior = 1, prior_type){
  s <- sample$var;
  adj <- sample$adj;
  p <- nrow(s);
  
  
  rowsum <- apply(sample$sample, 1, sum)
  data = sample$sample/rowsum
  output <- ccrepe(x = data, iterations = 1000, min.subj = 10)
  aa <- output$sim.score
  diag(aa) <-1
  pvalue <- output$p.values
  diag(pvalue) <- 1
  C_opt <- matrix(0, p, p);
  CI_opt <- matrix(0, p, p);
  
  
  R_opt <- aa;
  
  A_opt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(R_opt[i,j]) > 0 && pvalue[i,j] < 0.05  ){
        A_opt[i,j] <- 1;
        A_opt[j,i] <- 1;
        C_opt[i,j] <- R_opt[i,j];
        C_opt[j,i] <- R_opt[j,i];
        CI_opt[i,j] <- R_opt[i,j];
        CI_opt[j,i] <- R_opt[j,i];
      }
    }
  }
  diag(A_opt) <- 1;
  
  return(list("cor" = R_opt, "cov" = C_opt, "icov" = CI_opt, "adj" = A_opt ));
}

