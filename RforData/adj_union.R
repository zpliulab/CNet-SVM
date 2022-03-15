## INPUT gene and gene net£¬OUTPUT adj matrix

library(igraph)

rm(list=ls())

setwd('D:\\E\\...Data\\RTCGA')


# DATA ----------------------------------------------------------------------

gene <- read.csv('UNgene_component.csv')
genes <- data.frame(gene[,1])
net <- read.csv('UNgene_comp_net.csv')


G <- graph_from_data_frame(net, directed=F, vertices=genes)
print(G, e=TRUE, v=TRUE)



# Graph to adj matrix ---------------------------------------------------------------
# adj <- as_adjcaency_matrix(G,sparse=FALSE) 
adj <- get.adjacency(G,sparse=FALSE) 
View(adj[1:10,1:10])
# write.csv(adj, 'adjmatrix_comp_UNG.csv')


