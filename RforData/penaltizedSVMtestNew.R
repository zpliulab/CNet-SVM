

rm(list = ls())


library(caret)
library(dplyr)
## penalizedSVM 
library(mlegp)
library(e1071)
library(MASS)
library(tgp)  
library(maptree)
library(penalizedSVM)

library(pROC)

# load function -----------------------------------------------------------
source("D:\\E\\...\\RforData\\svmpenalized.R")

seed <- 123

# input data ----------------------------------------------------------------------
setwd('D:\\E\\...\\Data\\RTCGA')
x <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))


# method list --------------------------------------------------------------------
methodlist <- c("1norm", "DrHSVM", "scad", "scad+L2")

 
j <- 1

# "1norm", "DrHSVM(Elastic Net)", "scad", "scad+L2(Elastic SCAD)"-- svmfs
name <- methodlist[j]
result <- svmpenalized(data, name)

# save results ----------------------------------------------------------------------
featscad <- result[[1]]
aucscad <- result[[2]]
predscad <- result[[3]]

setwd('D:\\E\\...\\Data\\RTCGA\\result\\featurenew') 
save(featscad, file = paste("./RdataFeature/",paste(name,".RData")))
save(aucscad, file = paste("./RdataAUC/",paste(name,".RData")))
save(predscad, file = paste("./RdataPred/",paste(name,".RData")))



