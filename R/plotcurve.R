source('~/Desktop/HMP_project/my_R/evaluation.R');
load(file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/synthetic_new/scale_free_100_400.RData");

glasso_p_indicator = TRUE;
cclasso_indicator = TRUE;
REBACCA_indicator = TRUE;
SparCC_indicator = TRUE;
mb_indicator = TRUE;
glasso_indicator = TRUE;
ccrepe_indicator = TRUE;

if (glasso_p_indicator){
  templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  temppred <- split(t(predgp), rep(1:nrow(predgp), each = ncol(predgp)))
  glasso_p_ROC <- ROCnew(templabel, temppred)
  perfgp <- glasso_p_ROC$perf;
  aucgp <- glasso_p_ROC$auc;
  aucstdgp <- glasso_p_ROC$aucsd
  
  acc <- accEva(templabel, temppred);
  print("glasso_p:")
  print(sprintf("%.3f (%.3f)", mean(acc), sd(acc)));
  plot(perfgp,lwd=3,avg="threshold", col="red", yaxis.at=c(-1, 2.0), xaxis.at = c(-1, 2),ann = FALSE)
  axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), NA, cex.axis=.7, font=1, tck=.03)
  axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), NA, cex.axis=.7, font=1, tck=.03)
}
if (cclasso_indicator){
  templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  temppred <- split(t(predcc), rep(1:nrow(predcc), each = ncol(predcc)))
  cclasso_ROC <- ROCnew(templabel, temppred)
  perfcc <- cclasso_ROC$perf;
  auccc <- cclasso_ROC$auc;
  aucstdcc <- cclasso_ROC$aucsd
  acc <- accEva(templabel, temppred);
  print("cclasso:")
  print(sprintf("%.3f (%.3f)", mean(acc), sd(acc)));
  #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEcc), sd(RMSEcc), mean(RMSEFcc), sd(RMSEFcc), auccc, aucstdcc))
  #plot(perfcc,col="grey82",lty=3)
  plot(perfcc,lwd=3,avg="threshold",add=TRUE, col = "darkgreen", lty=2)
}

if (REBACCA_indicator){
  templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  temppred <- split(t(predr), rep(1:nrow(predr), each = ncol(predr)))
  REBACCA_ROC <- ROCnew(templabel, temppred)
  perfr <- REBACCA_ROC$perf;
  aucr <- REBACCA_ROC$auc;
  aucstdr <- REBACCA_ROC$aucsd
  acc <- accEva(templabel, temppred);
  print("REBACCA:")
  print(sprintf("%.3f (%.3f)", mean(acc), sd(acc)));
  #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEr), sd(RMSEr), mean(RMSEFr), sd(RMSEFr), aucr, aucstdr))
  #plot(perfr,col="grey82",lty=3)
  plot(perfr,lwd=3,avg="threshold",add=TRUE, col = "gold2", lty=2)
}

if (SparCC_indicator){
  templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  temppred <- split(t(preds), rep(1:nrow(preds), each = ncol(preds)))
  SparCC_ROC <- ROCnew(templabel, temppred)
  perfs <- SparCC_ROC$perf;
  aucs <- SparCC_ROC$auc;
  aucstds <- SparCC_ROC$aucsd
  acc <- accEva(templabel, temppred);
  print("SparCC:")
  print(sprintf("%.3f (%.3f)", mean(acc), sd(acc)));
  #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEs), sd(RMSEs), mean(RMSEFs), sd(RMSEFs), aucs, aucstds))
  #plot(perfs,col="grey82",lty=3)
  plot(perfs,lwd=3,avg="threshold",add=TRUE, col = "blue", lty=2)
}

if (mb_indicator){
  templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  temppred <- split(t(predm), rep(1:nrow(predm), each = ncol(predm)))
  mb_ROC <- ROCnew(templabel, temppred)
  perfmb <- mb_ROC$perf;
  aucmb <- mb_ROC$auc;
  aucstdmb <- mb_ROC$aucsd
  acc <- accEva(templabel, temppred);
  print("mb:")
  print(sprintf("%.3f (%.3f)", mean(acc), sd(acc)));
  #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(ICOVm), sd(ICOVm), mean(ICOVFm), sd(ICOVFm), aucmb, aucstdmb))
  #plot(perfmb,col="grey82",lty=3)
  plot(perfmb,lwd=3,avg="threshold",add=TRUE, col = "cyan3", lty=2)
}

if (glasso_indicator){
  templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  temppred <- split(t(predgl), rep(1:nrow(predgl), each = ncol(predgl)))
  gl_ROC <- ROCnew(templabel, temppred)
  perfgl <- gl_ROC$perf;
  aucgl <- gl_ROC$auc;
  aucstdgl <- gl_ROC$aucsd
  acc <- accEva(templabel, temppred);
  print("glasso:")
  print(sprintf("%.3f (%.3f)", mean(acc), sd(acc)));
  #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(ICOVgl), sd(ICOVgl), mean(ICOVFgl), sd(ICOVFgl), aucgl, aucstdgl))
  #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(RMSEgl), sd(RMSEgl), mean(RMSEFgl), sd(RMSEFgl), aucgl, aucstdgl))
  #plot(perfgl,col="grey82",lty=3)
  plot(perfgl,lwd=3,avg="threshold",add=TRUE, col = "gray8", lty=2)
}

if (ccrepe_indicator){
  templabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  temppred <- split(t(predre), rep(1:nrow(predre), each = ncol(predre)))
  re_ROC <- ROCnew(templabel, temppred)
  perfre <- re_ROC$perf;
  aucre <- re_ROC$auc;
  aucstdre <- re_ROC$aucsd
  acc <- accEva(templabel, temppred);
  print("ccrepe:")
  print(sprintf("%.3f (%.3f)", mean(acc), sd(acc)));
  #print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(ICOVre), sd(ICOVre), mean(ICOVFre), sd(ICOVFre), aucre, aucstdre))
  #plot(perfre,col="grey82",lty=3)
  plot(perfre,lwd=3,avg="threshold",add=TRUE, col = "purple", lty=2)
}
#save(sampleLabel, predcc, predgl, predgp, predm, predr, predre, preds, perfcc, perfgl, perfgp, perfmb, perfr, perfre, perfs, file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/synthetic/scale_free_100_200.RData")
