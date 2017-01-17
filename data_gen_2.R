OTUnum <- 50;

graph <- make_graph('cluster', OTUnum, OTUnum)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))