rm(list=ls())

library(pROC)
library(MASS)
library(glmnet)


working_dir = "D:\\E\\...NSLR"
setwd(working_dir)

source('Fun_Auxiliary.R')
source('Fun_NSLR.R')

# ------------------------------------------
# Generate a simple simulation data
snrlam0 = 3
f.num0 = 150
f.ture = 75
sim.data = GET.SIM.DATA2(smaple.num = 120, feature.num = f.num0, 
                         feature.ture = f.ture, random.seed = 260, 
                         snrlam=snrlam0)
adj = get_sim_prior_Net(f.num0, 80, 0.03, 0.01)
dim(adj)
# write.csv(adj, file = "adj_example.csv", row.names = F)


library(igraph)
g <- graph.adjacency(adj,mode="undirected")
plot(g)


# -------------------------------------------------------------------

Train.id = 1:84
Valid.id = 85:120
w.true = sim.data$w

# Training data
X1 = sim.data$X[Train.id,]; y1 =sim.data$y[Train.id]; 
train_data <- as.matrix(cbind(y1, X1))
colnames(train_data) <- c("Label", c(1:f.num0))
dim(train_data)
# View(train_data[,1:10])
write.table(train_data, file = "Data_train\\260.txt",quote = F, sep = "\t")
write.csv(train_data, file = "Data_train_matlab\\260.txt")


# Testing data
X2 = sim.data$X[Valid.id,]; y2 =sim.data$y[Valid.id]; 
test_data <- cbind(y2,X2)
colnames(test_data) <- c("Label", c(1:f.num0))
dim(test_data)
# View(test_data[,1])
write.table(test_data, file = "Data_test\\260.txt",quote = F, sep = "\t")
write.csv(test_data, file = "Data_test_matlab\\260.txt")

