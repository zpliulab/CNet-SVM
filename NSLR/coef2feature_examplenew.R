## Result of SVMmainexample using SVMexamplenew file on server, random seed=500
## throh = 1e-2


###############   theta_selection (Select features from the coefficient theta output by matlab)  ##############

rm(list=ls())

setwd('D:\\E\\...\\SVMexample\\Data_theta')                             


theta0 = read.csv('theta_500.csv', header=F, sep=",") [1,]
theta = as.matrix(read.csv('theta_500.csv', header=F, sep=",")[-1,])

gene <- as.matrix(c(1:150))
rownames(theta) <- gene[,1]


label <- which(abs(theta) >= 1e-3)
feature <- gene[label,1]


setwd('D:\\E\\...\\NSLR')
net_used <- as.matrix(read.csv("list.csv",header = T))
node_used <- as.matrix(feature)

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


nodename <- get.vertex.attribute(g)
feature_new <- as.matrix(nodename$name)
colnames(feature_new) <- c('biomarker')



clu <- components(p1)
groups(clu)

component <- groups(clu)$'1'
setwd('D:\\E\\...\\SVMexample\\Feature') 
# write.csv(component, "feature_500.csv", row.names = F)

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
setwd('D:\\E\\...\\SVMexample\\FeatureNet') 
# write.csv(ed,"featurenet_500.csv.csv",row.names = F,quote = F)
