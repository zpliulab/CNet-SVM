
library(mlegp)
library(e1071)
library(MASS)


library(tgp)  
library(maptree)
library(penalizedSVM)


rm(list = ls())

seed <- 123   

setwd('D:\\E\\...Data\\RTCGA') 



data_all <- t(as.matrix(read.table("TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt", header = T, check.names = FALSE, sep="\t")))
dim(data_all)    # 224 793
View(data_all[,1:10])
x_all <- data_all[,-1]    # sample*gene
dim(x_all)    # 224 792
class(x_all)
y_all <- data_all[,1]
y_all = 2*as.numeric(as.factor(y_all))-3 
dim(y_all)
class(y_all)


## train
data <- as.matrix(read.table("TCGA_pro_outcome_TN_log_comp_UNtrain.txt", header = T, check.names = FALSE, sep="\t"))
dim(data)    # 158 793
View(data[,1:10])
x <- data[,-1]    # sample*gene
class(x)
y <- data[,1]
y = 2*as.numeric(as.factor(y))-3 
dim(y)
class(y)


## test
data_test <- as.matrix(read.table("TCGA_pro_outcome_TN_log_comp_UNtest.txt", header = T, check.names = FALSE, sep="\t"))
x_test <- data_test[,-1]    # sample*gene
y_test <- data_test[,1]
y_test = 2*as.numeric(as.factor(y_test))-3 



# SAVE results --------------------------------------------------------------------
# "scad-SVM", "Lasso-SVM", "scad+L2(L2SCAD-SVM)", "DrHSVM(Enet-SVM)" -- svmfs

setwd('D:\\E\\²©Ê¿\\R_³ÌÐò\\SVM\\Data\\RTCGA\\result\\feature') 

# lambda, up and down, the same with paper ElasticSCAD 
bounds <- t(data.frame(log2lambda1 = c(-10, 10)))
colnames(bounds)<-c("lower", "upper")

# SCAD   -----------------------------------------------------------

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
View(feature_scad)
# write.csv(feature_scad, file = "feature_scad.csv", row.names = F)


# predict
test.error.scad <- predict(scad_fit, x_test, y_test)
print(test.error.scad$tab)
test.error.scad$sensitivity
test.error.scad$specificity

# lambda
print(paste("minimal 5-fold cv error -- SCAD method:", scad_fit$model$fit.info$fmin,
            "by log2(lambda1)=", scad_fit$model$fit.info$xmin))

print(" all lambdas with the same minimum? ")
print(scad_fit$model$fit.info$points.fmin)

print(paste(scad_fit$model$fit.info$neval, "visited points"))

print("overview: over all visitied points in tuning parameter space with corresponding cv errors")
print(data.frame(Xtrain = scad_fit$model$fit.info$Xtrain,
                 cv.error = scad_fit$model$fit.info$Ytrain))


.plot.EPSGO.parms (scad_fit$model$fit.info$Xtrain, 
                   scad_fit$model$fit.info$Ytrain,
                   bound = bounds, Ytrain.exclude=10^16, plot.name=NULL)



# lasso  -----------------------------------------------------------

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
# write.csv(feature_lasso, file = "feature_lasso.csv", row.names = F)

# predict
test.error.lasso <- predict(lasso_fit, x_test, y_test)
test.error.lasso <- predict(lasso_fit, x_all, y_all)
print(test.error.lasso$tab)
test.error.lasso$sensitivity
test.error.lasso$specificity

# lambda
print(paste("minimal 5-fold cv error -- Lasso method:", lasso_fit$model$fit.info$fmin,
            "by log2(lambda1)=", lasso_fit$model$fit.info$xmin))

print(" all lambdas with the same minimum? ")
print(lasso_fit$model$fit.info$points.fmin)

print(paste(lasso_fit$model$fit.info$neval, "visited points"))

print("overview: over all visitied points in tuning parameter space with corresponding cv errors")
print(data.frame(Xtrain = lasso_fit$model$fit.info$Xtrain,
                 cv.error = lasso_fit$model$fit.info$Ytrain))

.plot.EPSGO.parms (lasso_fit$model$fit.info$Xtrain, 
                   lasso_fit$model$fit.info$Ytrain,
                   bound = bounds, Ytrain.exclude=10^16, plot.name=NULL)


# Ealstic net -----------------------------------------------------------

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
# write.csv(feature_elastic, file = "feature_elastic.csv", row.names = F)

# predict
test.error.elastic <- predict(elastic_fit, x_test, y_test)
test.error.elastic <- predict(elastic_fit, x_all, y_all)
print(test.error.elastic$tab)
test.error.elastic$sensitivity
test.error.elastic$specificity

# lambda
print(paste("minimal 5-fold cv error -- Enet method:", elastic_fit$model$fit.info$fmin,
            "by log2(lambda1)=", elastic_fit$model$fit.info$xmin))

print(" all lambdas with the same minimum? ")
print(elastic_fit$model$fit.info$points.fmin)

print(paste(elastic_fit$model$fit.info$neval, "visited points"))

print("overview: over all visitied points in tuning parameter space with corresponding cv errors")
print(data.frame(Xtrain = elastic_fit$model$fit.info$Xtrain,
                 cv.error = elastic_fit$model$fit.info$Ytrain))

.plot.EPSGO.parms (elastic_fit$model$fit.info$Xtrain, 
                   elastic_fit$model$fit.info$Ytrain,
                   bound = bounds, Ytrain.exclude=10^16, plot.name=NULL)



# L2SCAD  -----------------------------------------------------------

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
# write.csv(feature_l2scad, file = "feature_l2scad.csv", row.names = F)

# predict
test.error.l2scad <- predict(l2scad_fit, x_test, y_test)
test.error.l2scad <- predict(l2scad_fit, x_all, y_all)
print(test.error.l2scad$tab)
test.error.l2scad$sensitivity
test.error.l2scad$specificity

# lambda
l2scad_fit$model$cv.info


