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
setwd('D:\\E\\²©Ê¿\\R_³ÌÐò\\SVM\\Data\\RTCGA')
x <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))
dim(data)    #  224 793
data[1:5,1:5]


solength <- 30
## first
set.seed(123*1)
training.samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
train.data  <- as.matrix(data[training.samples, ])
test.data <- as.matrix(data[-training.samples, ]) 

## train
x <- train.data[,-1]    # sample*gene
y <-  2*as.numeric(as.factor(train.data[,1]))-3 
x_test <- test.data[,-1]    # sample*gene
y_test = 2*as.numeric(as.factor(test.data[,1]))-3 


# SVM-RFE -----------------------------------------------------------------
# result <- svmRFE(train.data, k=1, halve.above=100)    # standard SVM-RFE, you can use k=1.
featureRankedList = svmrfeFeatureRanking(x,y)
rank <- featureRankedList[1:solength]
feature <- colnames(x)[rank]


# fit
svmfit = svm(x[ , featureRankedList[1:solength]], y, cost = 10, kernel="linear")  # linear svm
summary(svmfit)

# result
fea <- c()
fea <- c(fea, feature)

# predict
test <- predict(svmfit, x_test[ , featureRankedList[1:solength]])
pred <- cbind(y_test, as.numeric(test))
p <- plot.roc(y_test, as.numeric(test), print.auc=T)
auclist <- as.numeric(p$auc)



