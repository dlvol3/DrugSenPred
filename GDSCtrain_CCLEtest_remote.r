# DrugSenPred project   #### A shadow copy of the CCLE>=>GDSC version, this one,
# instead of using CCLE as training data, recruits GDSC to train the model.

## PCPV2 ready

# Author: Yue Zhang
# Create date: 03/Aug/2018
# Contact: yue.zhang@lih.lu


### Usage
# Input: CCLE: IC50, RPKM; GDSC: IC50, TPM
# Output: 
#       drugname.pdf: Training ROC, test ROC/PRC using the other dataset
#       Cellline_number.txt
#       SummaryofRandomForest for each drug in tab files
#       Result table of prediction for each drug in tab files
library(tidyverse)
library(h2o)
library(DMwR)
h2o.init(nthreads = -1)
############################################################
############################################################
### Stage 1 Data preparation(Expression data, RNA-Seq)######
############################################################
############################################################

#======================================================================================================
# CCLE RPKM>TPM re-formatting 
#======================================================================================================

# CCLE(RPKM)and GDSC(TPM)data tables
#getwd()
setwd('/mnt/pcpnfs/homedirs/yzhang/CCLEGDSC_RF/data/')
CCLE_raw <- read.delim("CCLE_RPKM.txt",header = TRUE,sep = '\t')
#CCLE_raw[1:4,1:4]


GDSC_raw <- read.delim("GDSC.txt",sep="\t",header=TRUE)


### Process the CCLE data

### gsub the version of the ENSG
CCLE_raw$Name <- gsub("\\..+","",CCLE_raw$Name)
CCLE_raw$Description <- gsub("\\..+","",CCLE_raw$Description)

names(CCLE_raw) <- gsub("^X","",colnames(CCLE_raw))
names(CCLE_raw) <- gsub("\\.\\..+","",colnames(CCLE_raw))

# Convert RPKM to TPM
#colSums(round(apply(CCLE_mid[,3:5],2,function(x) x/sum(x)*10e6),2) # test for converting RPKM to TPM

CCLE_exp <- cbind(CCLE_raw[,1:2],round(apply(CCLE_raw[,3:length(colnames(CCLE_raw))],2,function(x) x/sum(x)*10e6),1))
### CCLE TPM data well-prepared

#======================================================================================================
# GDSC data 
#======================================================================================================
## remove NAs
GDSC_raw[is.na(GDSC_raw)] <- 0

### gsub the name of the cell line  1. remove ..  2. gsub . to -

names(GDSC_raw) <- gsub("\\.\\..+","",colnames(GDSC_raw))

names(GDSC_raw) <- gsub("\\.","-",colnames(GDSC_raw))

names(GDSC_raw) <- gsub("^X","",colnames(GDSC_raw))   # remove the annoying X at beginning

## Remove non-overlap genes in GDSC(57503) based on genes in CCLE(56318)
GDSC_exp <- GDSC_raw[GDSC_raw$`Gene-ID` %in% CCLE_exp$Name,]  # An overlap of 51115 features has been found


## Remove non-overlap genes in CCLE(56318) based on genes in GDSC_exp(51115)
CCLE_exp <- CCLE_exp[CCLE_exp$Name %in% GDSC_exp$`Gene-ID`,] 


#======================================================================================================
# Data preparation stage 1 finished 
#======================================================================================================

############################################################
############################################################
##  Data preparation Stage 2, map to IC50,label classes  ###
############################################################
############################################################

library(data.table)

##transpose the CCLE expression data

CCLE_expt <- transpose(as.data.frame(CCLE_exp))   # Transpose
rownames(CCLE_expt) <- colnames(as.data.frame(CCLE_exp))
colnames(CCLE_expt) <- CCLE_exp$Name
CCLE_expR <- CCLE_expt[c(-1,-2),]   # remove the useless rows
CCLE_expR$Cellline_CCLE <- rownames(CCLE_expR)  # make a new column for left_join later


##transpose the GDSC expression data

GDSC_expt <- transpose(as.data.frame(GDSC_exp))   # Transpose
rownames(GDSC_expt) <- colnames(as.data.frame(GDSC_exp))
colnames(GDSC_expt) <- GDSC_exp$`Gene-ID`
GDSC_expR <- GDSC_expt[c(-1,-2),]   # remove the useless rows
GDSC_expR$Cellline_GDSC <- rownames(GDSC_expR)  # make a new column for left_join later

#======================================================================================================
# 'Transpose the exp data' finished! 
#======================================================================================================

#======================================================================================================
# Deal with the IC50 data 
#======================================================================================================

allIC50 <- read.delim("IC50_all_in_table.txt",header = TRUE, sep = '\t')

###########################################    IC50CCLE
IC50CCLE <- subset(allIC50, select = c("drug", "ccle.name", "CCLE_log10_IC50"))
apply(IC50CCLE, 2, function(x) any(is.na(x) | is.infinite(x)))  ## NA checker

## Remove cell line with NAs
IC50CCLE <- na.omit(IC50CCLE)
apply(IC50CCLE, 2, function(x) any(is.na(x) | is.infinite(x)))  ## NA checker


###########################################    IC50GDSC
IC50GDSC <- subset(allIC50, select = c("drug", "gdsc.name", "GDSC_log10_IC50"))
apply(IC50GDSC, 2, function(x) any(is.na(x) | is.infinite(x)))  ## NA checker

## Remove cell line with NAs
IC50GDSC <- na.omit(IC50GDSC)
apply(IC50GDSC, 2, function(x) any(is.na(x) | is.infinite(x)))  ## NA checker


#================================

####Check number of drugs, cell lines avaliable
CCLE_IC50_cellline_numbers <- IC50CCLE %>%
  group_by(drug) %>%
  summarise(n = n())

CCLE_IC50_drug_numbers <- IC50CCLE %>%
  group_by(ccle.name) %>%
  summarise(n = n())

GDSC_IC50_cellline_numbers <- IC50GDSC %>%
  group_by(drug) %>%
  summarise(n = n())

GDSC_IC50_drug_numbers <- IC50GDSC %>%
  group_by(gdsc.name) %>%
  summarise(n = n())

numberofavalibleDrugsCelllineList <- list(CCLE_IC50_cellline_numbers,CCLE_IC50_drug_numbers,
                                          GDSC_IC50_cellline_numbers,GDSC_IC50_drug_numbers)

cat(capture.output(print(numberofavalibleDrugsCelllineList),
                   file = paste0(getwd(),"/output/numberofavalibleDrugsCelllineList.txt")))

####============ Throw the terminal info into a .txt


#======================================================================================================
# Map the cutoff value back to the IC50 table to label the cellline-drug pairs 
#======================================================================================================

CCLEinflection <- read.delim("CCLE_inflection_point.txt", sep = '\t', header = TRUE)  # pre-prepared data table, from 'nature paper' S2

GDSCinflection <- read.delim("GDSC_inflection_point.txt", sep = '\t', header = TRUE)

### check the drug name
ot <- levels(as.factor(IC50CCLE$drug))
nt <- levels(as.factor(CCLEinflection$DRUG))
setdiff(ot,nt)   #####!!! Crizotinib in inflection data is actually the PF-2341066 in IC50 data..
## ref: https://en.wikipedia.org/wiki/Crizotinib

# step 1 remove the - in IC50 data
IC50CCLE$drug <- gsub("-","",IC50CCLE$drug)
IC50GDSC$drug <- gsub("-","",IC50GDSC$drug)

# step 2 remove the - in inflection data
CCLEinflection$DRUG <- gsub("-","",CCLEinflection$DRUG)
GDSCinflection$DRUG <- gsub("-","",GDSCinflection$DRUG)

# step 3 sub "PF2341066" into "Crizotinib" in the IC50 data
IC50CCLE$drug <- gsub("PF2341066","Crizotinib",IC50CCLE$drug)
IC50GDSC$drug <- gsub("PF2341066","Crizotinib",IC50GDSC$drug)

#### A Checker here again..
ot <- levels(as.factor(IC50CCLE$drug))
nt <- levels(as.factor(CCLEinflection$DRUG))
setdiff(ot,nt)    # Pass!
### CCLE first
IC50CCLE$SENRES <- ''
for (i in 1:length(rownames(IC50CCLE))){
  if (CCLEinflection[grepl(IC50CCLE$drug[i],CCLEinflection$DRUG),2] < IC50CCLE$CCLE_log10_IC50[i]){
    IC50CCLE$SENRES[i] <- "RES"
  }
  else if (CCLEinflection[grepl(IC50CCLE$drug[i],CCLEinflection$DRUG),2] >= IC50CCLE$CCLE_log10_IC50[i]){
    IC50CCLE$SENRES[i] <- "SEN"
  }
}

### then GDSC
IC50GDSC$SENRES <- ''
for (i in 1:length(rownames(IC50GDSC))){
  if (GDSCinflection[grepl(IC50GDSC$drug[i],GDSCinflection$DRUG),2] < IC50GDSC$GDSC_log10_IC50[i]){
    IC50GDSC$SENRES[i] <- "RES"
  }
  else if (GDSCinflection[grepl(IC50GDSC$drug[i],GDSCinflection$DRUG),2] >= IC50GDSC$GDSC_log10_IC50[i]){
    IC50GDSC$SENRES[i] <- "SEN"
  }
}

#======================================================================================================
# Labeling fishished 
#======================================================================================================
## Are the classes balanced or not?
# Checker

table(as.factor(IC50CCLE$SENRES))
table(as.factor(IC50GDSC$SENRES))

### Not very balanced...


####Explore the distribution of the IC50 values

# To be added

#ggplot(IC50CCLE, aes(y = IC50CCLE$CCLE_log10_IC50, x = IC50CCLE$drug))+
#  geom_violin(stat = "ydensity", position = "dodge")


####

############################################################
############################################################
######           looping for building the RF           #####
############################################################
############################################################

time.start <- proc.time()[3]
#### Before that..Let's use one drug to build RF for testing the scripts and the data prepraration

# Candidate: 17AAG   ####Test passed 06/Aug/2018

#======================================================================================================
# Loop of 15 drugs
#======================================================================================================
for (i in 1:length(levels(as.factor(IC50CCLE$drug)))){
  drugg <- levels(as.factor(IC50CCLE$drug))[i]
  
  pdf(paste0(getwd(),"/output/",drugg,"_GDSC.pdf"))  # The pdf receiving the output
  OnedrugIC <- subset(IC50CCLE, drug == drugg, select = c("ccle.name","SENRES"))
  
  OnedrugAll <- left_join(y= CCLE_expR, x = OnedrugIC, by = c( "ccle.name" = "Cellline_CCLE" ))
  
  
  table(apply(OnedrugAll, 1, function(x) any(is.na(x) | is.infinite(x))))  ## NA checker
  #
  ## 471 mapped, 31 are not#
  #
  # Remove NAs rows
  
  OnedrugAll <- na.omit(OnedrugAll)
  
  ## make sure there is no NAs
  table(apply(OnedrugAll, 1, function(x) any(is.na(x) | is.infinite(x))))  ## NA checker
  ## Pass!
  
  ### Remove the name of the cellline?
  OnedrugAll <- OnedrugAll[,-1]
  
  ### Make Chr into Num and Factor 06/08/2018
  OnedrugAll$SENRES <- as.factor(OnedrugAll$SENRES)
  OnedrugAll[,2:length(colnames(OnedrugAll))] <- sapply(OnedrugAll[,2:length(colnames(OnedrugAll))],as.double)
  
  ###sapply is much more faster than mutate_if!!!!!!
  
  #### SMOTE the data
  #OnedrugAllSMOTE <- SMOTE(SENRES ~.,OnedrugAll, perc.over = 2000, k = 5, perc.under = 100)
  #table(OnedrugAllSMOTE$SENRES)
  
  ### Check the columns data type
  #str(OnedrugAll[,1:5])
  
  #======================================================================================================
  # Prepare data for validation/test GDSC here
  #======================================================================================================
  OnedrugICtest <- subset(IC50GDSC, drug == drugg, select = c("gdsc.name","SENRES"))
  
  OneDrugTest <- left_join(y= GDSC_expR, x = OnedrugICtest, by = c( "gdsc.name" = "Cellline_GDSC" ))
  
  
  table(apply(OneDrugTest, 1, function(x) any(is.na(x) | is.infinite(x))))  ## NA checker
  #
  ## 471 mapped, 31 are not#
  #
  # Remove NAs rows
  
  OneDrugTest <- na.omit(OneDrugTest)
  
  ## make sure there is no NAs
  table(apply(OneDrugTest, 1, function(x) any(is.na(x) | is.infinite(x))))  ## NA checker
  ## Pass!
  
  ### Remove the name of the cellline?
  OneDrugTest <- OneDrugTest[,-1]
  
  ### Make Chr into Num and Factor 06/08/2018
  OneDrugTest$SENRES <- as.factor(OneDrugTest$SENRES)
  OneDrugTest[,2:length(colnames(OneDrugTest))] <- sapply(OneDrugTest[,2:length(colnames(OneDrugTest))],as.double)
  
  
  
  ########################
  # H2O.ai
  ########################
  
  set.seed(3)
  Onedrugh2ot <- as.h2o(OnedrugTest)
  OnedrugTestH2o <- as.h2o(OneDrugAll)
  
  # splits <- h2o.splitFrame(Onedrugh2ot,c(0.75),seed = 3)
  # train <- h2o.assign(splits[[1]],"train.hex")
  # test <- h2o.assign(splits[[2]],"test.hex")
  #str(train[,1:5])
  rf <- h2o.randomForest(
    training_frame = Onedrugh2ot,
    x = 2:51116,
    y = 1,
    ntrees = 1000,
    score_each_iteration = T,
    max_depth = 15,
    seed = 1234
  )                                 # Model training using CCLE
  
  
  
  #======================================================================================================
  # Result of the model/prediction/output the result datatable/draw the ROC and PRC 
  #======================================================================================================
  
  
  plot(h2o.performance(rf))    ### print the training ROC 
  #summary(rf)
  cat(capture.output(print(summary(rf)),
                     file = paste0(getwd(),"/output/summaryrfGDSC_",drugg,".txt")))   ### print summary into a txt
  
  
  result <- as.data.frame(h2o.predict(rf,OnedrugTestH2o))
  re <- as.data.frame(cbind(as.vector(OnedrugTestH2o$SENRES),result))
  colnames(re) <- c("Original","Predicted","Pred_RES","Pred_SEN")
  write.table(re,file = paste0(getwd(),"/output/PredictionGDSCtoCCLE_",drugg,".txt"),sep = '\t', col.names = TRUE,row.names = FALSE)
  ## write result table into a txt
  
  
  #======================================================================================================
  # ROC based on re
  #======================================================================================================
  
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
    ggtitle(paste0("ROC curve of", drugg, " AUC = ", auc))
  
  
  PRCc <- ggplot(plotready, aes(y = precision,x = TPR))+
    geom_path(color = "red", size = 2)+
    theme_bw()+
    ggtitle(paste0("PRC curve of ", drugg))
  
  
  print(ROCc)   ### print into pdf ROC
  print(PRCc)   ### print into pdf PRC
  
  ############################################################################
  
  dev.off()
  time.end <- time.start - proc.time()[3]
  
  message("*****",drugg,", Random forest training finished in ", round(time.end/60,2)," min.")
}


#### CCLE>=>=>GDSC finished





