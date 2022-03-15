
rm(list = ls())
library(igraph)

setwd('D:\\E\\...NSLR') 

num = 150
gene <- as.matrix(c(1:num))
genes <- data.frame(gene[,1])

net <- read.csv('list.csv')
g <- graph_from_data_frame(net, directed=F)



## breadth-first search
diameter(g) 

## TRUE, the diameters of the connected components
diameter(g, unconnected=TRUE)

## FALSE, the number of vertices
diameter(g, unconnected=FALSE)

## returns a path with the actual diameter
get_diameter(g) 

## returns two vertex ids, connected by the diameter path.
farthest_vertices(g) 
# diameter(g, weights=NA)

# cut  ----------------------------------------------------------------

## TRUE, only the minumum cut value is returned
min_cut(g, source = "132", target = "90", value.only = TRUE)


library(igraph)
g <- graph_from_data_frame(net, directed=F)

node <- as.matrix(get_diameter(g))
lab1 <- as.matrix(rownames(node))
colnames(lab1) <- c("node")


my_vector <- function(gene, lab1){
  k <- length(as.matrix(gene))
  vec1 <- matrix(0,1,k)
  
  l <- length(lab1)

  for (i in 1:l) {
    # i= 1
    vec1[which(gene[,1]==lab1[i])] <- c("-1")
  }
  
  x_first <- which(gene[,1]==lab1[1])
  x_end <- which(gene[,1]==lab1[l])
  vec1[,c(x_first,x_end)] <- c("1")
  
  return(vec1)
}

vec2 <- my_vector(gene, lab1)
vec3 <- t(vec2)
View(vec3)
# write.table(vec3,file = "vector_hat_example.txt", quote=F, sep="\t", row.names = F)

