
###############   theta_selection  ##############

rm(list=ls())

setwd('D:\\E\\...Matlab_code\\SVM\\Data_theta')                             

theta0 = read.csv('theta_lambda45.csv', header=F, sep=",") [1,]
theta = read.csv('theta_hat_lambda45.csv', header=F, sep=",")


setwd('D:\\...Data\\RTCGA')  
gene <- read.csv('UNgene_component.csv', header = T, sep=',')
rownames(theta) <- gene[,1]


label <- which(abs(theta) >= 1e-5)
feature <- gene[label,1]


# feature in net ------------------------------------------------------------

net_used <- as.matrix(read.csv("UNgene_comp_net.csv",header = T))
node_used <- as.matrix(feature)

k1 <- which(net_used[,1] %in% node_used) 
k2 <- which(net_used[,2] %in% node_used) 
length(intersect(k1,k2))   
used <- net_used[intersect(k1,k2),]


library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP)
ed <- as_edgelist(p1, names = TRUE)

g <- p1
plot(g, layout=layout.fruchterman.reingold, 
     vertex.size = 8,  
     vertex.label = V(g)$name, 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 0.1, 
     vertex.label.color = "black"  
)


# gene in Graoh as biomarker  ----------------------------------------------------
nodename <- get.vertex.attribute(g)
feature_new <- as.matrix(nodename$name)
colnames(feature_new) <- c('biomarker')


# max net
clu <- components(p1)
groups(clu)
component <- groups(clu)$'1'
# write.csv(component, file = "result\\feature\\feature_component.csv", row.names = F)

node_used <- as.matrix(component)
k1 <- which(net_used[,1] %in% node_used) 
k2 <- which(net_used[,2] %in% node_used) 
length(intersect(k1,k2))   
used <- net_used[intersect(k1,k2),]


library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP)
# p1 <- simplify(PP,remove.multiple = TRUE, remove.loops = TRUE) 
ed <- as_edgelist(p1, names = TRUE)

g <- p1
plot(g, layout=layout.fruchterman.reingold,  
     vertex.size = 8,   
     vertex.label = V(g)$name,  
     vertex.label.cex = 0.8,  
     vertex.label.dist = 0.1,  
     vertex.label.color = "black"   
)
# write.csv(ed,"result\\feature_net.csv",row.names = F,quote = F)


# biomarker label ----------------------------------------------------------

label <- which(rownames(theta) %in% component)


# predict ---------------------------------------------------------------

theta_hat <- as.matrix(theta[label,])
rownames(theta_hat) <- c(gene[label,1])
class(theta_hat)
theta_hat <- as.matrix(theta_hat)


library(mlegp)
library(e1071)
library(MASS)
library(penalizedSVM)
library(pROC)
require(caret)


setwd('D:\\E\\...Data\\RTCGA')  
data <- t(as.matrix(read.table('TCGA_pro_outcome_TN_log_comp_UNtest.txt', header=T, sep='\t', check.names = F)))
data <- t(data)


xtrain <- data[,-1]    # sample*gene
# change class labels to -1 and 1
ytrain = 2*as.numeric(as.factor(data[,1]))-3 


datanew <- data[,-1]
View(as.matrix(datanew)[,label, drop=FALSE][,1:10])
sep = as.matrix(datanew)[,label, drop=FALSE] %*% theta_hat + theta0
pred.class = factor(2*as.numeric (sep > 1) -1)
tab.classes <- table(ytrain, pred.class)
plot.roc(as.numeric(ytrain), as.numeric(sep), print.auc=T)
confusionMatrix(tab.classes, positive = "1") 


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

