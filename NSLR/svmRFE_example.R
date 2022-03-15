###################  SVM-RFE ÔÚexampleÊý¾Ý  #############################


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
setwd('D:\\E\\...\\NSLR')

solength <- 30
train.data  <- as.matrix(read.table("Data_train\\10.txt", header = T, check.names = FALSE))
test.data <- as.matrix(read.table("Data_test\\10.txt", header = T, check.names = FALSE)) 

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
setwd('D:\\E\\...\\Data\\ExampleData') 
write.csv(feature, file = "feature_svmrfe10.csv", row.names = F)

# fit
svmfit = svm(x[ , featureRankedList[1:solength]], y, cost = 10, kernel="linear")  # linear svm
summary(svmfit)


# predict
test <- predict(svmfit, x_test[ , featureRankedList[1:solength]])
pred <- cbind(y_test, as.numeric(test))
p <- plot.roc(y_test, as.numeric(test), print.auc=T)
auclist <- as.numeric(p$auc)


# index
pred.class = factor(2*as.numeric (test > 0) -1)
tab.classes <- table(pred.class, pred[,1])

tp <- tab.classes[2, 2]
tn <- tab.classes[1, 1]
fp <- tab.classes[2, 1]
fn <- tab.classes[1, 2]


accuracy <- (tp + tn)/(tp + tn + fp + fn)
precision <- (tp)/(tp+fp)
sensitivity <- tp/(tp + fn)
specificity <- tn/(tn + fp)
F_measure <- 2*precision*sensitivity/(precision+sensitivity)   

tab.classes
accuracy 
precision
sensitivity 
specificity
F_measure
