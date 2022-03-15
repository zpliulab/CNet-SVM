###################  SVM-RFE  #############################


rm(list = ls())

library(caret)
library(dplyr)
library(e1071)
library(pROC)

setwd('D:\\E\\...\\RforData') 
# source('msvmRFE.R')   # form  package  mSVM-RFE   2005
source('svmrfeFeatureRanking.R')    # form  package  sigFeature   2020
source('svmpenalizedsvmRFE.R')

# input data ----------------------------------------------------------------------
setwd('D:\\E\\...\\Data\\RTCGA')
x <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))

result <- svmpenalizedsvmRFE(data)

# save results ----------------------------------------------------------------------
featscad <- result[[1]]
aucscad <- result[[2]]
predscad <- result[[3]]


name <- "SVM-RFE"
setwd('D:\\E\\...\\Data\\RTCGA\\result\\featurenew') 
save(featscad, file = paste("./RdataFeature/",paste(name,".RData")))
save(aucscad, file = paste("./RdataAUC/",paste(name,".RData")))
save(predscad, file = paste("./RdataPred/",paste(name,".RData")))
