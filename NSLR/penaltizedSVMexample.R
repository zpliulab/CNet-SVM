## 2021.8.23 改进 2011_Elastic SCAD as a novel penalization method 的 Rpackage
## lambda序列的设置非常重要，目前l1-SVM的序列有问题
## 2021.9.19 改进。使用区间搜索寻优，替换原来的网格法，看一下迭代的次数10是否有问题

###################################################################################

## penalizedSVM 是依赖于下面3个packages的
# install.packages("mlegp")
library(mlegp)
library(e1071)
library(MASS)


library(tgp)  # 本地安装
# install.packages("maptree")
library(maptree)
library(penalizedSVM)


rm(list = ls())

seed <- 123    # svmfs 函数必须的

setwd('D:\\E\\博士\\R_程序\\SVM\\R\\NSLR') 


# 读入数据 --------------------------------------------------------------------

## 训练
data <- as.matrix(read.table("Data_train\\100.txt", header = T, check.names = FALSE, sep="\t"))
dim(data)    # 158 793
# View(data[,1:10])
x <- data[,-1]    # sample*gene
class(x)
y <- data[,1]
y = 2*as.numeric(as.factor(y))-3 
dim(y)
class(y)


## 测试
data_test <- as.matrix(read.table("Data_test\\100.txt", header = T, check.names = FALSE, sep="\t"))
x_test <- data_test[,-1]    # sample*gene
y_test <- data_test[,1]
y_test = 2*as.numeric(as.factor(y_test))-3 



# 保存结果 --------------------------------------------------------------------
# "scad", "1norm", "scad+L2(Elastic SCAD)", "DrHSVM(Elastic Net)" -- svmfs

setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\ExampleData') 

# lambda 的上下界，与论文中的一样
bounds <- t(data.frame(log2lambda1 = c(-10, 10)))
colnames(bounds)<-c("lower", "upper")

# SCAD 模型 -----------------------------------------------------------
 
# computation intensive
print("start interval search -- SCAD method")

scad_fit <- svmfs(x, y, fs.method="scad",
                  cross.outer = 0, 
                  bounds = bounds,
                  # grid.search = "discrete",
                  # lambda1.set =lambda,
                  # parms.coding = "none",
                  parms.coding = "log2",  
                  show = "none",
                  maxIter = 10, 
                  inner.val.method = "cv", 
                  cross.inner = 5,
                  maxevals = 500,
                  seed = seed, 
                  verbose = FALSE) 

system.time(scad_fit)
print("scad final model")
print(str(scad_fit$model))

# result
w_scad <- scad_fit$model$w
b_scad <- scad_fit$model$b
fea_num_scad <- scad_fit$model$xind
feature_scad <- colnames(x)[fea_num_scad]
# View(feature_scad)
# write.csv(feature_scad, file = "feature_scad100.csv", row.names = F)

# predict
test.error.scad <- predict(scad_fit, x_test, y_test)
# test.error.scad <- predict(scad_fit, x_all, y_all)
print(test.error.scad$tab)
tab.classes <- test.error.scad$tab 
test.error.scad$sensitivity
test.error.scad$specificity
library(pROC)
plot.roc(as.numeric(y_test), as.numeric(test.error.scad$fitted), print.auc=T)


tp <- tab.classes[2, 2]
tn <- tab.classes[1, 1]
fp <- tab.classes[2, 1]
fn <- tab.classes[1, 2]


accuracy <- (tp + tn)/(tp + tn + fp + fn)
precision <- (tp)/(tp+fp)
sensitivity <- tp/(tp + fn)
specificity <- tn/(tn + fp)
F_measure <- 2*precision*sensitivity/(precision+sensitivity)   


accuracy 
precision
sensitivity 
specificity
F_measure


# 最优参数
print(paste("minimal 5-fold cv error -- SCAD method:", scad_fit$model$fit.info$fmin,
            "by log2(lambda1)=", scad_fit$model$fit.info$xmin))

print(" all lambdas with the same minimum? ")
print(scad_fit$model$fit.info$points.fmin)

print(paste(scad_fit$model$fit.info$neval, "visited points"))

print("overview: over all visitied points in tuning parameter space with corresponding cv errors")
print(data.frame(Xtrain = scad_fit$model$fit.info$Xtrain,
                 cv.error = scad_fit$model$fit.info$Ytrain))

# create 3 plots on one screen:
# 1st plot: distribution of initial points in tuning parameter space
# 2nd plot: visited lambda points vs. cv errors
# 3rd plot: the same as the 2nd plot, Ytrain.exclude points are excluded.
# The value cv.error = 10^16 stays for the cv error for an empty model !
# pdf(file = "figure_scad.pdf",width = 10,height = 6)
.plot.EPSGO.parms (scad_fit$model$fit.info$Xtrain, 
                   scad_fit$model$fit.info$Ytrain,
                   bound = bounds, Ytrain.exclude=10^16, plot.name=NULL)
# dev.off()


# lasso 模型 -----------------------------------------------------------

# computation intensive
print("start interval search -- Lasso method")

lasso_fit <- svmfs(x, y, fs.method="1norm",
                  cross.outer = 0, 
                  bounds = bounds,
                  # grid.search = "discrete",
                  # lambda1.set =lambda,
                  # parms.coding = "none",
                  parms.coding = "log2",  
                  show = "none",
                  maxIter = 10, 
                  inner.val.method = "cv", 
                  cross.inner = 5,
                  maxevals = 500,
                  seed = seed, 
                  verbose = FALSE) 

system.time(lasso_fit)
print("lasso final model")
print(str(lasso_fit$model))

# result
w_lasso <- lasso_fit$model$w
b_lasso <- lasso_fit$model$b
fea_num_lasso <- lasso_fit$model$xind
feature_lasso <- colnames(x)[fea_num_lasso]
# write.csv(feature_lasso, file = "feature_lasso100.csv", row.names = F)

# predict
test.error.lasso <- predict(lasso_fit, x_test, y_test)
# test.error.lasso <- predict(lasso_fit, x_all, y_all)
print(test.error.lasso$tab)
tab.classes <- test.error.lasso$tab
test.error.lasso$sensitivity
test.error.lasso$specificity
plot.roc(as.numeric(y_test), as.numeric(test.error.lasso$fitted), print.auc=T)


tp <- tab.classes[2, 2]
tn <- tab.classes[1, 1]
fp <- tab.classes[2, 1]
fn <- tab.classes[1, 2]


accuracy <- (tp + tn)/(tp + tn + fp + fn)
precision <- (tp)/(tp+fp)
sensitivity <- tp/(tp + fn)
specificity <- tn/(tn + fp)
F_measure <- 2*precision*sensitivity/(precision+sensitivity)   


accuracy 
precision
sensitivity 
specificity
F_measure


# 最优参数
print(paste("minimal 5-fold cv error -- Lasso method:", lasso_fit$model$fit.info$fmin,
            "by log2(lambda1)=", lasso_fit$model$fit.info$xmin))

print(" all lambdas with the same minimum? ")
print(lasso_fit$model$fit.info$points.fmin)

print(paste(lasso_fit$model$fit.info$neval, "visited points"))

print("overview: over all visitied points in tuning parameter space with corresponding cv errors")
print(data.frame(Xtrain = lasso_fit$model$fit.info$Xtrain,
                 cv.error = lasso_fit$model$fit.info$Ytrain))

# create 3 plots on one screen:
# 1st plot: distribution of initial points in tuning parameter space
# 2nd plot: visited lambda points vs. cv errors
# 3rd plot: the same as the 2nd plot, Ytrain.exclude points are excluded.
# The value cv.error = 10^16 stays for the cv error for an empty model !
# pdf(file = "figure_lasso.pdf",width = 10,height = 6)
.plot.EPSGO.parms (lasso_fit$model$fit.info$Xtrain, 
                   lasso_fit$model$fit.info$Ytrain,
                   bound = bounds, Ytrain.exclude=10^16, plot.name=NULL)
# dev.off()


# Ealstic net 模型 -----------------------------------------------------------

# computation intensive
print("start interval search -- Enet method")

elastic_fit <- svmfs(x, y, fs.method="DrHSVM",
                   cross.outer = 0, 
                   bounds = bounds,
                   # grid.search = "discrete",
                   # lambda1.set =lambda,
                   # parms.coding = "none",
                   parms.coding = "log2",  
                   show = "none",
                   maxIter = 10, 
                   inner.val.method = "cv", 
                   cross.inner = 5,
                   maxevals = 500,
                   seed = seed, 
                   verbose = FALSE) 

system.time(elastic_fit)
print("elastic final model")
print(str(elastic_fit$model))

# result
w_elastic <- elastic_fit$model$w
b_elastic <- elastic_fit$model$b
fea_num_elastic <- elastic_fit$model$xind
feature_elastic <- colnames(x)[fea_num_elastic]
# write.csv(feature_elastic, file = "feature_elastic100.csv", row.names = F)

# predict
test.error.elastic <- predict(elastic_fit, x_test, y_test)
# test.error.elastic <- predict(elastic_fit, x_all, y_all)
print(test.error.elastic$tab)
tab.classes <- test.error.elastic$tab 
test.error.elastic$sensitivity
test.error.elastic$specificity
plot.roc(as.numeric(y_test), as.numeric(test.error.elastic$fitted), print.auc=T)


tp <- tab.classes[2, 2]
tn <- tab.classes[1, 1]
fp <- tab.classes[2, 1]
fn <- tab.classes[1, 2]


accuracy <- (tp + tn)/(tp + tn + fp + fn)
precision <- (tp)/(tp+fp)
sensitivity <- tp/(tp + fn)
specificity <- tn/(tn + fp)
F_measure <- 2*precision*sensitivity/(precision+sensitivity)   


accuracy 
precision
sensitivity 
specificity
F_measure


# 最优参数
print(paste("minimal 5-fold cv error -- Enet method:", elastic_fit$model$fit.info$fmin,
            "by log2(lambda1)=", elastic_fit$model$fit.info$xmin))

print(" all lambdas with the same minimum? ")
print(elastic_fit$model$fit.info$points.fmin)

print(paste(elastic_fit$model$fit.info$neval, "visited points"))

print("overview: over all visitied points in tuning parameter space with corresponding cv errors")
print(data.frame(Xtrain = elastic_fit$model$fit.info$Xtrain,
                 cv.error = elastic_fit$model$fit.info$Ytrain))

# create 3 plots on one screen:
# 1st plot: distribution of initial points in tuning parameter space
# 2nd plot: visited lambda points vs. cv errors
# 3rd plot: the same as the 2nd plot, Ytrain.exclude points are excluded.
# The value cv.error = 10^16 stays for the cv error for an empty model !
# pdf(file = "figure_elastic.pdf",width = 10,height = 6)
.plot.EPSGO.parms (elastic_fit$model$fit.info$Xtrain, 
                   elastic_fit$model$fit.info$Ytrain,
                   bound = bounds, Ytrain.exclude=10^16, plot.name=NULL)
# dev.off()




# Ealstic SCAD 模型 -----------------------------------------------------------

# computation intensive
print("start interval search -- L2SCAD method")


lambda <- c (seq(0.01 ,0.05, .01), seq(0.1,0.5, 0.2), 1 )
lambda <-lambda[2:3]

l2scad_fit <- svmfs(x, y, fs.method="scad+L2",
                    cross.outer = 0, 
                    grid.search = "discrete",
                    lambda1.set = lambda,
                    parms.coding = "none",
                    # bounds = bounds,
                    # parms.coding = "log2",
                    show = "none",
                    maxIter = 10, 
                    inner.val.method = "cv", 
                    cross.inner = 5,
                    maxevals = 500,
                    seed = seed, 
                    verbose = FALSE) 

system.time(l2scad_fit)
print("l2scad final model")
print(str(l2scad_fit$model))

# result
w_l2scad <- l2scad_fit$model$w
b_l2scad <- l2scad_fit$model$b
fea_num_l2scad <- l2scad_fit$model$xind
feature_l2scad <- colnames(x)[fea_num_l2scad]
# write.csv(feature_l2scad, file = "feature_l2scad100.csv", row.names = F)

# predict
test.error.l2scad <- predict(l2scad_fit, x_test, y_test)
# test.error.l2scad <- predict(l2scad_fit, x_all, y_all)
print(test.error.l2scad$tab)
tab.classes <- test.error.l2scad$tab
test.error.l2scad$sensitivity
test.error.l2scad$specificity
plot.roc(as.numeric(y_test), as.numeric(test.error.l2scad$fitted), print.auc=T)


tp <- tab.classes[2, 2]
tn <- tab.classes[1, 1]
fp <- tab.classes[2, 1]
fn <- tab.classes[1, 2]


accuracy <- (tp + tn)/(tp + tn + fp + fn)
precision <- (tp)/(tp+fp)
sensitivity <- tp/(tp + fn)
specificity <- tn/(tn + fp)
F_measure <- 2*precision*sensitivity/(precision+sensitivity)   


accuracy 
precision
sensitivity 
specificity
F_measure



# 最优参数
l2scad_fit$model$cv.info

# print(paste("minimal 5-fold cv error -- Elastic SCAD method:", l2scad_fit$model$fit.info$fmin,
#             "by log2(lambda1)=", l2scad_fit$model$fit.info$xmin))
# 
# print(" all lambdas with the same minimum? ")
# print(l2scad_fit$model$fit.info$points.fmin)
# 
# print(paste(l2scad_fit$model$fit.info$neval, "visited points"))
# 
# print("overview: over all visitied points in tuning parameter space with corresponding cv errors")
# print(data.frame(Xtrain = l2scad_fit$model$fit.info$Xtrain,
#                  cv.error = l2scad_fit$model$fit.info$Ytrain))
# 
# # create 3 plots on one screen:
# # 1st plot: distribution of initial points in tuning parameter space
# # 2nd plot: visited lambda points vs. cv errors
# # 3rd plot: the same as the 2nd plot, Ytrain.exclude points are excluded.
# # The value cv.error = 10^16 stays for the cv error for an empty model !
# # pdf(file = "figure_l2scad.pdf",width = 10,height = 6)
# .plot.EPSGO.parms (l2scad_fit$model$fit.info$Xtrain, 
#                    l2scad_fit$model$fit.info$Ytrain,
#                    bound = bounds, Ytrain.exclude=10^16, plot.name=NULL)
# # dev.off()
