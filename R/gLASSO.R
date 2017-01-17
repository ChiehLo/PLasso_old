library(glasso)
source('~/Desktop/HMP_project/SpiecEasi-master/R/community_graph.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/fitdistr.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/mvdistributions.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/normalization.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/plotNet.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/roc.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/spaRcc.R')
source('~/Desktop/HMP_project/SpiecEasi-master/R/SparseICov.R')



# set.seed(100)
# x<-matrix(rnorm(50*20),ncol=20)
# s<- var(x)
rho1 <- matrix(0.01, nrow = 50, ncol = 50)
rho1[3, 4] <- 1;
rho1[4, 3] <- 1;
#a<-glasso(s, rho=.01)
# a<-glasso(s, rho=rho1)
# icov <- a$wi;

library(R.matlab)
rawData = readMat("/Users/chiehlo/Desktop/HMP_project/syn_8.mat");
x <- rawData$Sample;
R <- rawData$R;
# P <- rawData$P;
A <- rawData$A;
data <- t(clr(x+1, 1))
s <- var(data)
rho2 <- rawData$rho;
a<-glasso(s, rho=rho2)
C_opt <- a$w;


OTUnum <- nrow(C_opt);
d <- 0;
for(i in 1:OTUnum){
  d[i] <- 1/sqrt(C_opt[i,i]);
}
D <- diag(d);
R_opt <- D%*%C_opt%*%D;
Error <- abs(R - R_opt);
RMSE <- sum(sum(Error))/(OTUnum*(OTUnum-1));

## glasso
#rho2 <- rawData$rho;
a<-glasso(s, rho=0.01)
C_opt <- a$w;


OTUnum <- nrow(C_opt);
d <- 0;
for(i in 1:OTUnum){
  d[i] <- 1/sqrt(C_opt[i,i]);
}
D <- diag(d);
R_opt <- D%*%C_opt%*%D;
Error <- abs(R - R_opt);
RMSE_g <- sum(sum(Error))/(OTUnum*(OTUnum-1));


# aa<-glasso(s,rho=.02, w.init=a$w, wi.init=a$wi)
# # example with structural zeros and no regularization,
# s=c(10,1,5,4,10,2,6,10,3,10)
# S=matrix(0,nrow=4,ncol=4)
# S[row(S)>=col(S)]=s
# S=(S+t(S))
# diag(S)<-10
# zero<-matrix(c(1,3,2,4),ncol=2,byrow=TRUE)
# a<-glasso(S,0,zero=zero)