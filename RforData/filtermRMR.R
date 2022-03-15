############################## filter mthod ###########################################


rm(list = ls())

# install.packages("mRMRe")
library(mRMRe)
library(caret)
library(dplyr)
library(pROC)


# load function -----------------------------------------------------------
source("D:\\E\\...\\RforData\\svmpenalizedmRMR.R")

# input data ----------------------------------------------------------------------
setwd('D:\\E\\...\\Data\\RTCGA')
x <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt", header = T, check.names = FALSE)
dim(x)
data <- data.frame(t(x))   
result <- svmpenalizedmRMR(data)


# save results ----------------------------------------------------------------------
featscad <- result[[1]]
aucscad <- result[[2]]
predscad <- result[[3]]


name <- "mRMR-SVM"
setwd('D:\\E\\...\\Data\\RTCGA\\result\\featurenew') 
save(featscad, file = paste("./RdataFeature/",paste(name,".RData")))
save(aucscad, file = paste("./RdataAUC/",paste(name,".RData")))
save(predscad, file = paste("./RdataPred/",paste(name,".RData")))









