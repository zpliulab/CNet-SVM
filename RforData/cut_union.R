# INPUT gene and gene net£¬OUTPUT cut node and cut_vector_UNgene


rm(list = ls())

library(igraph)


setwd('D:\\E\\...Data\\RTCGA')
gene <- read.csv('UNgene_component.csv')
genes <- data.frame(gene[,1])


net <- read.csv('UNgene_comp_net.csv')
g <- graph_from_data_frame(net, directed=F)



## dim£¬breadth-first search
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

# minimum cut  ----------------------------------------------------------------

## FALSE, the edges in the cut and a the two (or more) partitions are also returned.
min_cut(g, source = "ZBTB33", target = "OVCH2", value.only = FALSE)

## TRUE, only the minumum cut value is returned
min_cut(g, source = "ZBTB33", target = "OVCH2", value.only = TRUE)


dg <- as.directed(g)
min_cut(dg, value.only = FALSE)


st_min_cuts(dg, source = "ZBTB33", target = "OVCH2")
cut <- st_min_cuts(dg, source = "ZBTB33", target = "OVCH2")


cut$value
E(dg)[cut$cuts[[1]]]
V(dg)[cut$partition1s[[1]]]
cut$cuts[[2]] 
cut$partition1s[[2]]


# Minimum size vertex separators ------------------------------------------
min_separators(g)


# matrix and vector -----------------------------------------------------------------
setwd('D:\\E\\...Data\\RTCGA')


library(igraph)
gene <- read.csv('UNgene_component.csv',header = T)
net <- read.csv('UNgene_comp_net.csv')
g <- graph_from_data_frame(net, directed=F)


node <- as.matrix(get_diameter(g))
lab1 <- as.matrix(rownames(node))
colnames(lab1) <- c("node")


my_vector <- function(gene, lab1){
  k <- length(as.matrix(gene))
  vec1 <- matrix(0,1,k)
  l <- length(lab1)
  for (i in 2:l-1) {
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
# write.table(vec3,file = "cut_vector_UNgene.txt", quote=F, sep="\t", row.names = F)
