
rm(list=ls())
setwd('D:\\E\\...NSLR')
net <- as.matrix(read.csv("list.csv",header = T))
net[1,]

setwd('D:\\E\\...Data\\ExampleData')
node <- as.matrix(read.csv("feature_elastic500.csv",header = T))

node[1]
node_used <- node
net_used <- net
k1 <- which(net_used[,1] %in% node_used)   # 4701    2042
k2 <- which(net_used[,2] %in% node_used)   # 7133     957
length(intersect(k1,k2))    # 230  13      # 193
used <- net_used[intersect(k1,k2),]    # 192

library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP,remove.loops = T,remove.multiple = T)  # 最终的数对
ed <- as_edgelist(p1, names = TRUE)


clu <- components(p1)
groups(clu)

g <- p1


eb <- cluster_edge_betweenness(g)
eb[[1]]
eb[[2]]
eb[[3]]

n.vertices <- vcount(g); 
n.edges <- ecount(g);  

n.vertices
n.edges
