

# Feature selection -------------------------------------------------------



rm(list = ls())

# install.packages("SGL")

library(SGL)


setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\RTCGA')
x <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt", header = T, check.names = FALSE)
dim(x)
data <- data.frame(t(x))
# View(data[,1:10])
p <- dim(data)[2]-1


y <- as.matrix(t(x)[,-1])
dim(y)
y0 <- as.numeric(data[,1])

p <- dim(y)[2]


# #SGL
# data0 <- list(x=y,y=y0)
# 
# 
# set.seed(123)
# sgl <- cvSGL(data0 ,index <- c(rep(1:p,each=1)),min.frac=0.05)
# View(sgl$fit$beta)
# 
# c2 <- c(1:20)
# i <- 1
# 
# for(i in 1:20)
# {c2[i] <- length(colnames(y)[which(sgl$fit$beta[,i]!=0)])}
# 
# c2
# colnames(y[,which(sgl$fit$beta[,which(c2==71)]!=0)])




# machine learning --------------------------------------------------------

library(dplyr)
library(ROCR)
library(caret)
library(kernlab)


# class <- y0
# sglG <- colnames(y[,which(sgl$fit$beta[,which(c2==71)]!=0)])#12为SGL选择的基因个数

# svmpenalizedmRMR <- function(data){
condset <- c(3,13,14,17:20,23:25,28,32,33,36,38:41,44,46)

## first
set.seed(123*20)
training_samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
##  # sample*gene(col 1 -label)
train_data  <- as.matrix(data[training_samples, -1])  
test_data <- as.matrix(data[-training_samples, -1])  
## data.fram  --  train and test
train_data <- data.frame(train_data) 
test_data <- data.frame(test_data) 
y_train <- as.numeric(data[training_samples, 1])
y_test <- as.numeric(data[-training_samples, 1])





#SGL
data0 <- list(x=train_data,y=y_train)

sgl <- cvSGL(data0 ,index <- c(rep(1:p,each=1)),min.frac=0.05)
View(sgl$fit$beta)

c2 <- c(1:20)
for(i in 1:20)
{c2[i] <- length(colnames(y)[which(sgl$fit$beta[,i]!=0)])}

c2
sglG <- colnames(y[,which(sgl$fit$beta[,which(c2==28)]!=0)])
fea <- as.matrix(sglG)
colnames(fea) <- "SGLasso-SVM"

# write.csv(fea, file="result\\featurenew\\CsvdataFeatureSeleOnce\\SGL.csv",row.names = F)



# SGL result -------------------------------------------------------------

train  <- data.frame(data[training_samples, sglG])  
test <- data.frame(data[-training_samples, sglG])  


A1 <- ksvm(y_train~.,data=train, kernel="rbfdot",prob.model=TRUE)
result <-predict(A1,test,type="response")
list<-cbind(result,test$class)
# pred<-prediction(predictions=result[,2],labels=test$class)
pred<-prediction(predictions=result[,1],labels=y_test)
perd<-performance(pred,measure="tpr",x.measure="fpr")
plot(perd,main="ROC curve for SMS spam filter",col="blue",lwd=2)
abline(a=0,b=1,lwd=2,lty=2)
perf.auc<-performance(pred,measure="auc")
#str(perf.auc)
#unlist(perf.auc@y.values)
auc <- unlist(perf.auc@y.values)
table(result,y_test)
agreement <- result==y_test
f1 <- table(agreement)






