
rm(list = ls())

# install.packages("mRMRe")
library(mRMRe)
library(caret)
library(pROC)

seed <- 123    # svmfs need

# input data ----------------------------------------------------------------------
setwd('D:\\E\\...\\NSLR')
train.data  <- as.matrix(read.table("Data_train\\10.txt", header = T, check.names = FALSE))
test.data <- as.matrix(read.table("Data_test\\10.txt", header = T, check.names = FALSE)) 

## data.fram  --  train and test
train_data <- data.frame(train.data) 
test_data <- data.frame(test.data) 

mrmre_train <- mRMR.data(data=train_data[,-1], strata = factor(train_data[,1]))

## train -- classical mRMR feature selection
solength <- 30
mdata <- mrmre_train
classic.m <- mRMR.classic(mdata, feature_count=solength, target_indices=1)
## unable to find an inherited method for function ¡®featureNames¡¯ for signature 
sefs_rmre <- mdata@feature_names[solutions(classic.m)[[1]]]

# fit ---------------------------------------------------------------------
ens1_fs <- as.matrix(sefs_rmre)
## ens1_fs 
model <- apply(ens1_fs, 2, function(x, y) {
  ff <- as.formula(sprintf("%s ~ %s", colnames(y)[1], paste(x, collapse=" + ")))
  mm <- lm(formula=ff, data=y, model=FALSE)
  return(mm)
}, y=train_data)


# mRMR result -------------------------------------------------------------

feature <- sefs_rmre
setwd('D:\\E\\...\\Data\\ExampleData') 
write.csv(feature, file = "feature_mRMRSVM10.csv", row.names = F)

# predict ---------------------------------------------------------------------
test <- predict(object=model[[1]], newdata=test_data, type="response")
pred <- cbind(test_data[,1], as.numeric(test))
p <- plot.roc(test_data[,1], as.numeric(test), print.auc=T)

# index
pred.class = ifelse(test > 0.5, 1, 0)
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
