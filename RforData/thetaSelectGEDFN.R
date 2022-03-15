
## use GCWs method get top 1% genes, repeat 10 times, make Union

rm(list=ls())

setwd('D:\\E\\...Data\\RTCGA')
gene = read.csv(file = "allgene_component.csv", header=T, sep=",")

setwd('D:\\E\\...GEDFN-master')
myfile = list.files("rank")                
dir = paste("./rank/", myfile, sep="")     
n = length(dir)                            
theta <- as.matrix(read.table(file = dir[1], header=F))
rownames(theta) <- gene[,1]


ave <- as.matrix(quantile(theta[,1], seq(0.01, 1, 0.01)))[99]
label <- which(theta[,1] >= ave)  
feature <- gene[label,1]


theta = read.csv(file = "rank\\var_impo_svm1.csv",header=F, sep=",")  

for (i in 2:n){
  new_theta = read.csv(file = dir[i], header=F, sep=",")
  theta = cbind(theta,new_theta)

  ave <- as.matrix(quantile(theta[,i], seq(0.01, 1, 0.01)))[99]
  label <- which(theta[,i] >= ave) 
  feature <- cbind(feature, gene[label,1])
  }                                                

## union gene
feature1 <- feature[,1]
for (j in 2:dim(feature)[2]) {
  # j <- 4
  feature1 <- union(feature1, feature[,j]) 
}
unionfeature <-  feature1


setwd('D:\\E\\...Data')
# write.csv(unionfeature,file = "GEDFN169.csv", row.names = F)

