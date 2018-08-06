############################################################
############################################################
######            ROC/PRC based on re         ##############
############################################################
############################################################

# Author: Yue Zhang
# Create date: 05/Aug/2018
# Contact: yue.zhang@lih.lu

reAUC <- re[,c(-2,-4)]

nnumber <- 500
TPR <- rep(0,nnumber)
FPR <- rep(0,nnumber)
precision <- rep(0,nnumber)
cutpoint <- rep(0,nnumber)

p <- which(reAUC$Original == "RES")
n <- which(reAUC$Original == "SEN")
### Calculate the TPR/FPR

for(i in 1:nnumber){
  cutpoint[i] <- as.numeric(i/nnumber)
  predictionR <- reAUC$Pred_RES > cutpoint[i]
  pp <- intersect(which(predictionR == TRUE),p)   ### TP
  fp <- intersect(which(predictionR == TRUE),n)  ### Number of FP
  TPR[i] <- length(pp)/length(p)
  FPR[i] <- length(fp)/length(n)
  precision[i] <- length(pp)/length(which(predictionR == TRUE))
}

### Calculate the AUC

pos.value <- reAUC$Pred_RES[which(reAUC$Original == "RES")]
neg.value <- reAUC$Pred_RES[which(reAUC$Original == "SEN")]

aucs1 <- replicate(2000,mean(sample(pos.value,1000,replace = T) > sample(neg.value,1000,replace = T)))
auc <- round(mean(aucs1),4)

plotready <- cbind(FPR,TPR,precision)
colnames(plotready) <- c("FPR","TPR","precision")

plotready <- as.data.frame(plotready)

ROCc <- ggplot(plotready, aes(x = FPR,y = TPR))+
  geom_path(color = "red", size = 2)+
  geom_abline(intercept = 0, slope = 1, color='grey',size = 0.7)+
  theme_bw()+
  ggtitle(paste0("ROC curve, AUC = ", auc))


PRCc <- ggplot(plotready, aes(y = precision,x = TPR))+
  geom_path(color = "red", size = 2)+
  theme_bw()+
  ggtitle("PRC curve")





