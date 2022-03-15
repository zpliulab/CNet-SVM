## clear
rm(list = ls())

## package
library(caret)
library(dplyr)

# dir ----------------------------------------------------------------------
setwd('D:\\E\\...Data\\RTCGA')
dir.create("Data_train")
dir.create("Data_test")

# load data ----------------------------------------------------------------------
x <- read.table("TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))
# View(data[,1:10])

# data split -----------------------------------------------------------

## set 20 seed
condset <- c(3,13,14,17:20,23:25,28,32,33,36,38:41,44,46)

## first
i <- 1
set.seed(123*i)
training_samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
train_data  <- data[training_samples, ] 
test_data <- data[-training_samples, ] 

name <- as.character(i)
train <- as.matrix(train_data)
test <- as.matrix(test_data)

path <- paste("./Data_train/",paste(name,".txt",sep=""),sep="")
write.table(train,path,quote = F)

path <- paste("./Data_test/",paste(name,".txt",sep=""),sep="")
write.table(test,path,quote = F)
## second
for (i in condset) {
  print(paste("******* i= ********",i))
  set.seed(123*i)
  training_samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
  train_data  <- as.matrix(data[training_samples, ])  
  test_data <- as.matrix(data[-training_samples, ])  
  
  name <- as.character(i)
  train <- as.matrix(train_data)
  test <- as.matrix(test_data)
  
  path <- paste("./Data_train/",paste(name,".txt",sep=""),sep="")
  write.table(train,path,quote = F)
  
  path <- paste("./Data_test/",paste(name,".txt",sep=""),sep="")
  write.table(test,path,quote = F)
}
