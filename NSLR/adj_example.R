
rm(list=ls())
library(data.table)

# DATA --------------------------------------------------------------------

setwd('D:\\E\\...NSLR')   
adj <- as.matrix(read.csv('adj_example.csv'))

num = 150
label <- c(1:num)

colnames(adj) <- label
rownames(adj) <- label

adj[upper.tri(adj)] <- 0
adj_up <- adj


mat_df = function(m){
  stopifnot(is.matrix(m))
  res = as.data.frame.table(m)
  setnames(res,old = names(res),new = c("row","col","value"))
  res
}



tdf_list = mat_df(adj_up)
dim(tdf_list)   
head(tdf_list)

list1 <- tdf_list[-which(tdf_list[,3] == "0"),]
list2 <- list1[,-3]
# write.csv(list2, file = "list.csv", row.names = F, quote = F)
