
rm(list = ls())


# load Rdata --------------------------------------------------------------
setwd('D:\\E\\...\\Data\\RTCGA\\result\\featurenew') 
myfile = list.files("RdataFeature")    
dir = paste("./RdataFeature/", myfile, sep="")    
n = length(dir) 


## dir
dir
i <- 2
dir[i]


methodlist <- c("1norm", "CNet-SVM", "DrHSVM", "mRMR-SVM", "scad", "scadL2", "SVM-RFE")
name <- methodlist[i]
load(file = dir[i])


## count
numall <- matrix(nrow = 1,ncol=21)
for (i in 1:n) {
  num <- NULL
  for (j in 1:21) {
    name <- methodlist[i]
    load(file = dir[i])
    feature <- as.matrix(unlist(featscad[j]))
    colnames(feature) <- name
    p <-length(feature)
    num <- c(num, p)
  }
  numall <- rbind(numall, t(as.matrix(num)))
}
View(as.matrix(num))
View(numall)
nummat <- numall[-1,] 
rownames(nummat) <- methodlist
colnames(nummat) <- rep(1:21,1)


## Select the result of the 8th experiment
for (i in 1:n) {
  name <- methodlist[i]
  load(file = dir[i])
  feature <- as.matrix(unlist(featscad[8]))
  colnames(feature) <- name
  write.csv(feature, file = paste("./CsvdataFeatureSeleOnce/",paste(name,".csv")), row.names = F)
}

## end end  















## frequency of features
features <- featscad[c(1:20)]
# Sorted selection frequency across all 30 repetitions:
sort(table(unlist(features)), decreasing = TRUE)
feature <- as.data.frame(sort(table(unlist(features)), decreasing = TRUE))
colnames(feature) <- c(name, "Fren")
# write.csv(feature, file = paste("./CsvdataFeature/",paste(name,".csv")), row.names = F)


## stability 
library(stabm)
library(ggdendro)

stability <- NULL
for (i in 1:n) {
  load(file = dir[i])
  p <- length(unique(unlist(featscad[1:(length(featscad)-1)], use.names = FALSE)))
  stab <- stabilityHamming(featscad[1:(length(featscad)-1)], p)
  stability <- c(stability, stab)
  stability <- as.matrix(stability)
  colnames(stability) <- c("Stability")
}
rownames(stability) <- methodlist
stab0 <- t(stability)
# write.csv(stab0, file = "stability.csv", row.names = F)


## Average number of features
numave  <- NULL
for (k in 1:n) {
  # k = 2
  load(file = dir[k])
  features <- featscad[c(1:20)]
  num <- NULL
  for (j in 1:20) {
    num <- c(num, length(features[[j]]))
  }
  numave <- c(numave, apply(as.matrix(num), 2, mean))
}

numave <- as.matrix(numave)
rownames(numave) <- methodlist
numavet <- t(numave)
# write.csv(numavet, file = "NumResult.csv", row.names = F)




rm(list = ls())
# load Rdata --------------------------------------------------------------
setwd('D:\\E\\...\\RTCGA\\result\\featurenew') 
myfile = list.files("CsvdataFeature")    
dir = paste("./CsvdataFeature/", myfile, sep="")     
n = length(dir) 

# Who has the least features -------------------------------------------------------------
feanum <- NULL
for (i in 1:n) {
  num <- dim(as.matrix(read.csv(file = dir[i], header=T, sep=",")))[1]
  feanum <- c(feanum, num)
}
p <- min(feanum)

# output large list -------------------------------------------------------------------
fea <- as.matrix(read.csv(file = dir[1], header=T, sep=",")[1:p,])

for (i in 2:n){
  feanew = read.csv(file = dir[i], header=T, sep=",")[1:p,]
  fea = cbind(fea, feanew)
} 
fea
# write.csv(fea, file = "FsResult.csv", row.names = F)

## Selected in each method
dim(fea)
Reduce(intersect, fea[c(1,3,5,7,9,11,13)] )
Reduce(intersect, fea[c(3,1)] )
Reduce(intersect, fea[c(3,5)] )
Reduce(intersect, fea[c(3,7)] )
Reduce(intersect, fea[c(3,9)] )
Reduce(intersect, fea[c(3,11)] )
Reduce(intersect, fea[c(3,13)] )

method7 <- fea[c(1,3,5,7,9,11,13)]
# write.csv(method7, "methodfeature.csv", row.names = F)


# AUC ---------------------------------------------------------------------

rm(list = ls())


## load Rdata  
setwd('D:\\E\\²©Ê¿\\...\\Data\\RTCGA\\result\\featurenew') 
myfile = list.files("RdataAUC")    
dir = paste("./RdataAUC/", myfile, sep="")      
n = length(dir) 

## dir
dir
i <- 2
dir[i]

methodlist <- c("1norm", "CNet-SVM","DrHSVM", "mRMR-SVM", "scad", "scadL2", "SVM-RFE")
name <- methodlist[i]
load(file = dir[i])

View(aucscad)

load(file = dir[1])
value <- apply(as.matrix(aucscad), 2, mean)
for (i in 2:n){
  load(file = dir[i])
  valuenew <- apply(as.matrix(aucscad), 2, mean)
  value = cbind(value, valuenew)
}
value
colnames(value) <- methodlist
# write.csv(value, file = "AUCResult.csv", row.names = F)









