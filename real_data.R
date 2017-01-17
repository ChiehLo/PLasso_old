readOTU <- function(level = "genus", count_path, association = NULL, interaction = NULL){
  if(level == "genus"){
    rawData <- read.csv(count_path, sep = "\t", header=FALSE)
    rawData <- as.matrix(rawData)
    #rawData = readMat("/Users/chiehlo/Desktop/HMP_project/my_R/real_data_stool_genus.mat");
    sample = count_data_process(t(rawData), association);
    #sample = rawData;
  }
  else if(level == "species"){
    rawData <- t(as.matrix(read.table(count_path)));
    sample = fraction_data_process(rawData, association, interaction);
  }
  return(sample);
}

count_data_process <- function(rawData, association){
  p <- ncol(rawData);
  data <- t(clr(rawData+0.1, 1));
  s <- var(data);
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(s[i,i]);
  }
  D <- diag(d);
  R <- D%*%s%*%D;
  A <- association;
  return(list("sample" = rawData, "fraction" = data , "var" = s, "cor" = R, "adj" = association));
}

fraction_data_process <- function(rawData, association, interaction){
  p <- ncol(rawData);
  data <- t(clr(rawData+0.001, 1));
  s <- var(data);
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(s[i,i]);
  }
  D <- diag(d);
  R <- D%*%s%*%D;
  A <- association;
  return(list("sample" = rawData, "fraction" = data , "var" = s, "cor" = R, "adj" = association, "inter" = interaction));
}
