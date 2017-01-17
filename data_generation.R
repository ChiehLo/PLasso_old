## source library 

library(MASS)
library(Matrix)
library(matrixcalc)
library(glasso)
library(ROCR)
source('~/Desktop/HMP_project/SpiecEasi-master/R/community_graph.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/fitdistr.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/mvdistributions.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/normalization.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/plotNet.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/roc.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/spaRcc.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/SparseICov.R')

OTUnum <- 50; # number of OTU
Simulations <- 1e2; # number of sample
lambda <- 0.08; # glasso_p
lambda2 = 0.08; # glasso
prior <- 1;
numIt <- 1;
p <- OTUnum;
A <- matrix(0,p,p); # True Adjacency matrix
rho <- lambda*matrix(1, p, p); #gLasso regularizer
min_eigenvalue <- 1e-6; # make the covariance matrix PD
Sig <- matrix(0, p, p); # correlation matrix
Strength <- 0.15;
for (i in 1:p){
  for (j in i:p){
    r <- runif(1, min=0, max=1);
    if (r > 0.7){
      A[i,j] <- 1;
      A[j,i] <- 1;
      s <- runif(1, min=0, max=1);
      if (s > 0.5){
        Sig[i,j] = Strength;
        Sig[j,i] = Strength;
      }
      else{
        Sig[i,j] = -Strength;
        Sig[j,i] = -Strength;
      }
    }
    if (i==j){
      A[i,j] <- 1;
    }

  }
}

# regularization
for (i in 1:p){
  for (j in 1:p){
    if (abs(Sig[i,j]) == 0 && runif(1, min=0, max=1) > prior){
      rho[i,j] <- 10;
      rho[j,i] <- 10;
    }
  }
}
# Set the diagonal elements for covariance so that it’s positive definite
diag(Sig) <- max(1, min_eigenvalue - min(0, eigen(Sig)$values));
# Rescale covariance so that the diagonals are 1s
d2_Sig <- 1/sqrt(diag(Sig));
Sig2  <-  (d2_Sig)*Sig*(d2_Sig);
print(is.positive.definite(Sig2));
if(is.positive.definite(Sig2) == FALSE){
  Sigtemp <- nearPD(Sig2, keepDiag = TRUE)$mat;
  Sig3 <- matrix(Sigtemp@x, Sigtemp@Dim);
}else if(is.positive.definite(Sig2) == TRUE){
  Sig3 <- Sig2;
}

tempX = matrix(0, numIt+2,1);
tempY = matrix(0, numIt+2,1);
tempX[1] = 0;
tempY[1] = 0;
tempX[numIt+2] = 1;
tempY[numIt+2] = 1;
#generate sample
for (k in 1:numIt) {
  
Mu <- runif(p, 0, 1) - 0.5;
Sample <- exp(mvrnorm( n = Simulations, Mu , Sig3 ));
Sample <- round(Sample);
x <- Sample; # count data
R <- Sig3; #since we set the variance to 1
data <- t(clr(x+0.1, 1)) # log ratio transform
s <- var(data)

# sparcc
sparcc.amgut <- sparcc(x) #need to put count data
sparcc_R <- sparcc.amgut$Cor
Error_SparCC <- abs(R - sparcc_R);
RMSE_SpacCC <- sum(sum(Error_SparCC))/(OTUnum*(OTUnum-1));

#glasso_p
rho2 <- rho;
a<-glasso(s, rho=rho2)
C_opt <- a$w;
CI_opt <- a$wi;
d <- 0;
for(i in 1:OTUnum){
  d[i] <- 1/sqrt(C_opt[i,i]);
}
D <- diag(d);
R_opt <- D%*%C_opt%*%D;
Error <- abs(R - R_opt);
RMSE_Lasso_P <- sum(sum(Error))/(OTUnum*(OTUnum-1));

quan <- quantile(abs(R_opt), c(.25, .50, .80));

max_R <- max(max(R_opt));
A_opt <- matrix(0, nrow = OTUnum, ncol = OTUnum)
# for(i in 1:OTUnum){
#   for (j in i:OTUnum){
#     if (abs(R_opt[i,j]) > 0.01){
#       A_opt[i,j] <- 1;
#       A_opt[j,i] <- 1;
#     }
#   }
# }
for(i in 1:OTUnum){
  for (j in i:OTUnum){
    if (abs(CI_opt[i,j]) > 0){
      A_opt[i,j] <- 1;
      A_opt[j,i] <- 1;
    }
  }
}
print((nnzero(A)-p)/2);
ErEdge <- sum(sum(abs(A - A_opt)));
print(ErEdge);
count_pos <- 0;
count_fn <- 0;
count_fp <- 0;
for (i in 1:p){
  for (j in i:p){
    if (A[i,j] == 1 && A_opt[i,j] ==1 && i!=j){
      count_pos <- count_pos + 1;
    }
    if (A[i,j] == 1 && A_opt[i,j] == 0 && i!=j){
      count_fn <- count_fn + 1;
    }
    if (A[i,j] == 0 && A_opt[i,j] == 1 && i!=j){
      count_fp <- count_fp + 1;
    }
  }
}
A_pred <- A_opt[upper.tri(A_opt,diag=FALSE)];
A_label <- A[upper.tri(A,diag=FALSE)];
pred <- prediction( A_pred, A_label);
perf <- performance(pred, "tpr", "fpr");
auc <- performance(pred, 'auc')@y.values[[1]];
tempX[k+1] <- perf@x.values[[1]][2];
tempY[k+1] <- perf@y.values[[1]][2];
tempY <- tempY[order(tempX, decreasing=FALSE)]
tempX <- sort(tempX, decreasing=F);
#plot(perf);
}
perf@x.values[[1]] = tempX;
perf@y.values[[1]] = tempY;
plot(perf);

## gLasso 
# rho3 <- lambda2*matrix(1, p, p);
# a<-glasso(s, rho=rho3)
# C_opt <- a$w;
# d <- 0;
# for(i in 1:OTUnum){
#   d[i] <- 1/sqrt(C_opt[i,i]);
# }
# D <- diag(d);
# R_opt <- D%*%C_opt%*%D;
# Error <- abs(R - R_opt);
# RMSE_Lasso <- sum(sum(Error))/(OTUnum*(OTUnum-1));
# 
# max_R <- max(max(R_opt));
# A_opt <- matrix(0, nrow = OTUnum, ncol = OTUnum)
# for(i in 1:OTUnum){
#   for (j in i:OTUnum){
#     if (abs(R_opt[i,j]) > 1e-2){
#       A_opt[i,j] <- 1;
#       A_opt[j,i] <- 1;
#     }
#   }
# }
# ErEdge <- sum(sum(abs(A - A_opt)));
# count_pos_g <- 0;
# count_fn_g <- 0;
# count_fp_g <- 0;
# for (i in 1:p){
#   for (j in i:p){
#     if (A[i,j] == 1 && A_opt[i,j] ==1 && i!=j){
#       count_pos_g <- count_pos_g + 1;
#     }
#     if (A[i,j] == 1 && A_opt[i,j] == 0 && i!=j){
#       count_fn_g <- count_fn_g + 1;
#     }
#     if (A[i,j] == 0 && A_opt[i,j] == 1 && i!=j){
#       count_fp_g <- count_fp_g + 1;
#     }
#   }
# }
# # evaluation
# A_pred_g <- A_opt[upper.tri(A_opt,diag=FALSE)];
# A_label <- A[upper.tri(A,diag=FALSE)];
# pred_g <- prediction( A_pred_g, A_label);
# perf_g <- performance(pred_g, "tpr", "fpr");
# auc_g <- performance(pred_g, 'auc')@y.values[[1]];
#plot(perf_g);

