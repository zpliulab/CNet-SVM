

# Feature selection -------------------------------------------------------



rm(list = ls())

# install.packages("SGL")

library(SGL)


setwd('D:\\E\\博士\\R_程序\\SVM\\R\\NSLR')

# solength <- 30
train.data  <- as.matrix(read.table("Data_train\\500.txt", header = T, check.names = FALSE))
test.data <- as.matrix(read.table("Data_test\\500.txt", header = T, check.names = FALSE)) 


y <- as.matrix(train.data[,-1])
dim(y)
p <- dim(y)[2]


train_data <- as.matrix(train.data[,-1])    # sample*gene
y_train <-  as.numeric(train.data[,1])
test_data <- as.matrix(test.data[,-1])    # sample*gene
y_test = as.numeric(test.data[,1])



## data.fram  --  train and test
train_data <- data.frame(train_data) 
# View(train_data[,1:10])
test_data <- data.frame(test_data) 


# machine learning --------------------------------------------------------

library(dplyr)
library(ROCR)
library(caret)
library(kernlab)


## first
set.seed(123)

#SGL
data0 <- list(x=train_data,y=y_train)

sgl <- cvSGL(data0 ,index <- c(rep(1:p,each=1)),min.frac=0.05)
# View(sgl$fit$beta)

c2 <- c(1:20)
for(i in 1:20)
{c2[i] <- length(colnames(y)[which(sgl$fit$beta[,i]!=0)])}

c2
sglG <- colnames(y[,which(sgl$fit$beta[,which(c2==72)]!=0)])
fea <- as.matrix(sglG)
colnames(fea) <- "SGLasso-SVM"

# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\ExampleData') 
# write.csv(fea, file = "feature_SGLSVM500.csv", row.names = F)



# SGL result -------------------------------------------------------------

library(stringr)
sele <- str_c("X", sglG)



train  <- data.frame(train_data[, sele])  
test <- data.frame(test_data[, sele])  


A1 <- ksvm(y_train~.,data=train, kernel="rbfdot",prob.model=TRUE)
result <-predict(A1,test,type="response")
list <- cbind(result,y_test)
# pred<-prediction(predictions=result[,2],labels=test$class)
pred<-prediction(predictions=result[,1],labels=y_test)
perd<-performance(pred,measure="tpr",x.measure="fpr")
plot(perd,main="ROC curve for SGLasso-SVM",col="blue",lwd=2)
abline(a=0,b=1,lwd=2,lty=2)
# x11();plot(perd,main="ROC curve for SGLasso-SVM",col="blue",lwd=2);abline(a=0,b=1,lwd=2,lty=2)
perf.auc<-performance(pred,measure="auc")
#str(perf.auc)
#unlist(perf.auc@y.values)
auc <- unlist(perf.auc@y.values)
table(result,y_test)
agreement <- result==y_test
f1 <- table(agreement)




# prediction --------------------------------------------------------------

# predict
library(pROC)
pre <- predict(A1,test,type="response")
pred <- cbind(y_test, as.numeric(pre))
p <- plot.roc(y_test, as.numeric(pre), print.auc=T)
auclist <- as.numeric(p$auc)


# index
pred.class = factor(as.numeric(pre > 0.5))
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
auc
