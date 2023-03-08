
## 2023.3.5 Add the NNet method, compare with CNet-SVM method

rm(list = ls())

## Package
library(caret)
library(lattice)
library(ggplot2)


## pathway
# filepath <- '/Users/lilingyu/E/PhD/'
filepath <- '/home/lly/'

# load data ---------------------------------------------------------------

setwd(paste(filepath,'R/SVM/Data/RTCGA',sep = ""))
data <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt", header = T, 
                   sep = "\t", check.names = FALSE)
eset <- data[-1,]


Control <- rfeControl(functions = caretFuncs, method = "cv",
                      verbose = FALSE , returnResamp = "final")

trControl1 <- trainControl( method = "cv",
                            classProbs=TRUE,
                            summaryFunction = twoClassSummary)

## BRCA -- 0; Normal -- 1
disease <- as.factor(c(rep("BRCA",112), rep("Normal",112)))
# data[1,]


# NN-RFE ------------------------------------------------------------------

rf2 <- rfe(t(eset), disease, sizes = c(886),
           rfeControl = Control, trControl = trControl1, method = "nnet",
           tuneGrid = expand.grid(size = c(8), decay = c(0.1)),
           maxit = 30, MaxNWts = 100000
)
feature_sele <-rf2$optVariables
write.table(feature_sele, file = "result/featurenew/CsvdataFeatureSeleOnce/ranklist_NNETrfe.txt", quote=F, sep="\t")
