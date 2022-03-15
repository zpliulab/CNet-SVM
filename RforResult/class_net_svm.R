rm(list=ls())

library(pROC)
## l2s  --  las -- ela -- sca -- com


setwd('D:\\E\\...\\Data\\RTCGA\\result\\featuredatanew')


# Data  = read.table("Feature10780_CNetSVM.txt ", header = T, check.names = FALSE)
# Data2 = read.table("Feature_tcga_10780_CNetSVM.txt", header = T, check.names = FALSE)

# Data  = read.table("Feature42568_CNetSVM.txt ", header = T, check.names = FALSE)
# Data2 = read.table("Feature_tcga_42568_CNetSVM.txt", header = T, check.names = FALSE)

# Data  = read.table("Feature10797_CNetSVM.txt ", header = T, check.names = FALSE)
# Data2 = read.table("Feature_tcga_10797_CNetSVM.txt", header = T, check.names = FALSE)

# Data  = read.table("Feature45827_CNetSVM.txt ", header = T, check.names = FALSE)
# Data2 = read.table("Feature_tcga_45827_CNetSVM.txt", header = T, check.names = FALSE)

# Data  = read.table("Feature65194_CNetSVM.txt ", header = T, check.names = FALSE)
# Data2 = read.table("Feature_tcga_65194_CNetSVM.txt", header = T, check.names = FALSE)

# Data  = read.table("Feature38959_CNetSVM.txt ", header = T, check.names = FALSE)
# Data2 = read.table("Feature_tcga_38959_CNetSVM.txt", header = T, check.names = FALSE)

# Data  = read.table("Feature21422_CNetSVM.txt ", header = T, check.names = FALSE)
# Data2 = read.table("Feature_tcga_21422_CNetSVM.txt", header = T, check.names = FALSE)

# Data  = read.table("Feature26910_CNetSVM.txt ", header = T, check.names = FALSE)
# Data2 = read.table("Feature_tcga_26910_CNetSVM.txt", header = T, check.names = FALSE)

Data  = read.table("Feature20437_CNetSVM.txt ", header = T, check.names = FALSE)
Data2 = read.table("Feature_tcga_20437_CNetSVM.txt", header = T, check.names = FALSE)



#####################################  svm model -- built model makes prediction #########
x.train <- data.frame(t(Data2)[,-1])
y.train <- t(Data2)[,1]
x.test <- data.frame(t(Data)[,-1])
y.test <- t(Data)[,1]

library(e1071)
set.seed(666) 
tuned <- tune.svm(x.train,y.train, gamma = 10^(-6:-1), cost = 10^(1:2)) # tune
summary (tuned) # to select best gamma and cost

model <- svm(x.train, y.train, kernel = "radial", cost = tuned$best.parameters$cost, gamma=tuned$best.parameters$gamma,  scale = FALSE)
summary(model)




# Predict -----------------------------------------------------------------
p_test <- predict(model, x.test, type = "response")
## ROC ÇúÏß
library(pROC)
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, y.test)
setwd('D:\\E\\...\\Data\\RTCGA\\result\\A_test')
write.csv(A_test,"A_test_20437.csv", row.names=F)
names(A_test)<- c("p", "outcome")
# jpeg(file = "pAUC_train.jpg")
plot.roc(A_test$outcome, A_test$p, print.auc=T)
# dev.off()

