occurance <- function(occu_path, total_paper = 615268){
  a <- read.csv(occu_path, sep = "\t", header=FALSE)
  CountTable <- as.matrix(a)
  n = nrow(CountTable);
  P_value = matrix(0, n, n);
  total_paper = sum(diag(CountTable))
  # temp <- CountTable;
  # diag(temp) <- 0;
  # m <- sum(sum(temp));
  print(total_paper)
  for (i in 1:n){
    for (j in 1:n){
      if (j>i){
        contingency <- matrix(c(CountTable[i,j], CountTable[j,j] , CountTable[i, i], total_paper), nrow = 2)
        P_value[i,j] <- fisher.test(contingency, alternative = "two.sided")$p.value;
      }
    }
  }
  
  No_association <- matrix(0,n,n);
  adjust <- round(p.adjust(P_value, method = "bonferroni"), 5)
  corrected <- matrix(adjust, ncol=n, nrow = n);
  for (i in 1:n){
    for (j in 1:n){
      if (j>i && corrected[i,j] == 1){
        No_association[i,j] <- 1;
        No_association[j,i] <- 1;
      }
    }
  }
  return(No_association)
}

interaction <- function(interaction_path){
  a <- read.csv(interaction_path, sep = "\t", header=FALSE)
  interaction <- as.matrix(a)
  return(interaction)
}
