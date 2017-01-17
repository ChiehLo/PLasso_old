graph_select <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15, type = "random"){
  if(type=="random"){
    sample <- random_model(OTUnum, Simulations, Strength);
    print("random");
  }
  else if(type == "AR4"){
    sample <- AR4_model(OTUnum, Simulations, Strength);
    print("AR4");
  }
  else if(type == "hub"){
    sample <- hub_model(OTUnum, Simulations, Strength);
    print("hub");
  }
  else if(type == "cluster"){
    sample <- cluster_model(OTUnum, Simulations, Strength);
    print("cluster");
  }
  else if(type=="scale-free"){
    sample <- scale_free_model(OTUnum, Simulations, Strength);
    print("scale-free");
  }
  return(sample);
}


AR4_model <- function(OTUnum = 50, Simulations = 1e2, MaxStrength = 0.15 ){

  # p <- OTUnum;
  # A <- matrix(0,p,p); # True Adjacency matrix
  # min_eigenvalue <- 1e-6; # make the covariance matrix PD
  # Sig <- matrix(0, p, p); # correlation matrix
  # for (i in 1:p){
  #   for (j in 1:p){
  #     if (abs(i-j) == 0){
  #       A[i,j] = 1;
  #       Sig[i,j] = 1;
  #     }
  #     else if(abs(i-j) == 1){
  #       A[i,j] = 1;
  #       Sig[i,j] = 0.4;
  #     }
  #     else if(abs(i-j)==2){
  #       A[i,j] = 1;
  #       Sig[i,j] = 0.2;
  #     }
  #     else if(abs(i-j)==3){
  #       A[i,j] = 1;
  #       Sig[i,j] = 0.2;
  #     }
  #     else if(abs(i-j) == 4){
  #       A[i,j] = 1;
  #       Sig[i,j] = 0.1;
  #     }
  #   }
  # }
  # # Set the diagonal elements for covariance so that it’s positive definite
  # diag(Sig) <- max(1, min_eigenvalue - min(0, eigen(Sig)$values));
  # d2_Sig <- 1/sqrt(diag(Sig)); # Rescale covariance so that the diagonals are 1s
  # Sig2  <-  (d2_Sig)*Sig*(d2_Sig);
  # print(is.positive.definite(Sig2));
  # if(is.positive.definite(Sig2) == FALSE){
  #   Sigtemp <- nearPD(Sig2, keepDiag = TRUE)$mat;
  #   Sig3 <- matrix(Sigtemp@x, Sigtemp@Dim);
  # }else if(is.positive.definite(Sig2) == TRUE){
  #   Sig3 <- Sig2;
  # }
  # #generate sample
  # Mu <- runif(p, 0, 1) - 0.5;
  # Sample <- exp(mvrnorm( n = Simulations, Mu , Sig3 ));
  # Sample <- round(Sample);
  # x <- Sample; # count data
  # R <- Sig3; #since we set the variance to 1
  # data <- t(clr(x+0.1, 1)) # log ratio transform
  # s <- var(data)
  # return(list("sample" = x, "fraction" = data, "var" = s, "cor" = R, "adj" = A));
  
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "band", v = 0.3, u = 0 , g= 4, verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

hub_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.2 ){
  
  # p <- OTUnum;
  # A <- matrix(0,p,p); # True Adjacency matrix
  # min_eigenvalue <- 1e-6; # make the covariance matrix PD
  # Sig <- matrix(0, p, p); # correlation matrix
  # HubNumber <- 3;
  # Hubs <- sample.int(p, HubNumber);
  # for (i in 1:HubNumber){
  #   for (j in i:p){
  #     if ( runif(1, 0, 1) > 0.3 && Hubs[i] != j ){
  #       A[Hubs[i], j] = 1;
  #       A[j, Hubs[i]] = 1;
  #       Sig[Hubs[i], j] = Strength;
  #       Sig[j, Hubs[i]] = Strength;
  #     }
  #   }
  # }
  # for (i in 1:HubNumber){
  #   for (j in i:p){
  #     if ( runif(1, 0, 1) > 0.8 && length(which(Hubs==i)) == 0){
  #       A[i, j] = 1;
  #       A[j, i] = 1;
  #       Sig[i, j] = Strength;
  #       Sig[j, i] = Strength;
  #     }
  #     else if(i==j){
  #       A[i,j] = 1;
  #     }
  #   }
  # }
  # # Set the diagonal elements for covariance so that it’s positive definite
  # diag(Sig) <- max(1, min_eigenvalue - min(0, eigen(Sig)$values));
  # d2_Sig <- 1/sqrt(diag(Sig)); # Rescale covariance so that the diagonals are 1s
  # Sig2  <-  (d2_Sig)*Sig*(d2_Sig);
  # print(is.positive.definite(Sig2));
  # if(is.positive.definite(Sig2) == FALSE){
  #   Sigtemp <- nearPD(Sig2, keepDiag = TRUE)$mat;
  #   Sig3 <- matrix(Sigtemp@x, Sigtemp@Dim);
  # }else if(is.positive.definite(Sig2) == TRUE){
  #   Sig3 <- Sig2;
  # }
  # #generate sample
  # Mu <- runif(p, 0, 1) - 0.5;
  # Sample <- exp(mvrnorm( n = Simulations, Mu , Sig3 ));
  # Sample <- round(Sample);
  # x <- Sample; # count data
  # R <- Sig3; #since we set the variance to 1
  # data <- t(clr(x+0.1, 1)) # log ratio transform
  # s <- var(data)
  # return(list("sample" = x, "fraction" = data, "var" = s, "cor" = R, "adj" = A));
  
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "hub", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
  
}

cluster_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  
  # p <- OTUnum;
  # A <- matrix(0,p,p); # True Adjacency matrix
  # min_eigenvalue <- 1e-6; # make the covariance matrix PD
  # Sig <- matrix(0, p, p); # correlation matrix
  # cluster <- 2;
  # node_per_cluster <- p/cluster;
  # permutation <- sample.int(p, p);
  # for (i in 1:cluster){
  #   for (j in 1:node_per_cluster){
  #     for (k in 1:node_per_cluster){
  #       r <- runif(1, min=0, max=1);
  #       if(r > 0.7){
  #         A[permutation[(i-1)*node_per_cluster+j], permutation[(i-1)*node_per_cluster+k]] <- 1;
  #         A[permutation[(i-1)*node_per_cluster+k], permutation[(i-1)*node_per_cluster+j]] <- 1;
  #         Sig[permutation[(i-1)*node_per_cluster+j], permutation[(i-1)*node_per_cluster+k]] <- Strength;
  #         Sig[permutation[(i-1)*node_per_cluster+k], permutation[(i-1)*node_per_cluster+j]] <- Strength;
  #       }
  #     }
  #     if (j==k){
  #       A[permutation[(i-1)*node_per_cluster+j],permutation[(i-1)*node_per_cluster+j]] <- 1;
  #       Sig[permutation[(i-1)*node_per_cluster+j],permutation[(i-1)*node_per_cluster+j]] <- Strength;
  #     }
  #   }
  # }
  # 
  # # Set the diagonal elements for covariance so that it’s positive definite
  # diag(Sig) <- max(1, min_eigenvalue - min(0, eigen(Sig)$values));
  # d2_Sig <- 1/sqrt(diag(Sig)); # Rescale covariance so that the diagonals are 1s
  # Sig2  <-  (d2_Sig)*Sig*(d2_Sig);
  # print(is.positive.definite(Sig2));
  # if(is.positive.definite(Sig2) == FALSE){
  #   Sigtemp <- nearPD(Sig2, keepDiag = TRUE)$mat;
  #   Sig3 <- matrix(Sigtemp@x, Sigtemp@Dim);
  # }else if(is.positive.definite(Sig2) == TRUE){
  #   Sig3 <- Sig2;
  # }
  # #generate sample
  # Mu <- runif(p, 0, 1) - 0.5;
  # Sample <- exp(mvrnorm( n = Simulations, Mu , Sig3 ));
  # Sample <- round(Sample);
  # x <- Sample; # count data
  # R <- Sig3; #since we set the variance to 1
  # data <- t(clr(x+0.1, 1)) # log ratio transform
  # s <- var(data)
  # return(list("sample" = x, "fraction" = data, "var" = s, "cor" = R, "adj" = A));
  
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "cluster", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

random_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  
  # p <- OTUnum;
  # A <- matrix(0,p,p); # True Adjacency matrix
  # min_eigenvalue <- 1e-6; # make the covariance matrix PD
  # Sig <- matrix(0, p, p); # correlation matrix
  # for (i in 1:p){
  #   for (j in i:p){
  #     r <- runif(1, min=0, max=1);
  #     if (r > 0.7){
  #       A[i,j] <- 1;
  #       A[j,i] <- 1;
  #       s <- runif(1, min=0, max=1);
  #       if (s > 0.5){
  #         Sig[i,j] = Strength;
  #         Sig[j,i] = Strength;
  #       }
  #       else{
  #         Sig[i,j] = -Strength;
  #         Sig[j,i] = -Strength;
  #       }
  #     }
  #     if (i==j){
  #       A[i,j] <- 1;
  #     }
  #     
  #   }
  # }
  # 
  # # Set the diagonal elements for covariance so that it’s positive definite
  # diag(Sig) <- max(1, min_eigenvalue - min(0, eigen(Sig)$values));
  # d2_Sig <- 1/sqrt(diag(Sig)); # Rescale covariance so that the diagonals are 1s
  # Sig2  <-  (d2_Sig)*Sig*(d2_Sig);
  # print(is.positive.definite(Sig2));
  # if(is.positive.definite(Sig2) == FALSE){
  #   Sigtemp <- nearPD(Sig2, keepDiag = TRUE)$mat;
  #   Sig3 <- matrix(Sigtemp@x, Sigtemp@Dim);
  # }else if(is.positive.definite(Sig2) == TRUE){
  #   Sig3 <- Sig2;
  # }
  # #generate sample
  # Mu <- runif(p, 0, 1) - 0.5;
  # Sample <- exp(mvrnorm( n = Simulations, Mu , Sig3 ));
  # Sample <- round(Sample);
  # x <- Sample; # count data
  # R <- Sig3; #since we set the variance to 1
  # data <- t(clr(x+0.1, 1)) # log ratio transform
  # s <- var(data)
  # return(list("sample" = x, "fraction" = data, "var" = s, "cor" = R, "adj" = A));
  
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "random", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

scale_free_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "scale-free", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}


