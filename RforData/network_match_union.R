## INPUT gene    OUTPUB pathway gene-gene

rm(list=ls())
setwd('D:\\E\\...Data')
net <- as.matrix(read.csv("Regnetwork_hum.csv",header = T))
net[1,]


setwd('D:\...Data\\RTCGA')
node <- as.matrix(read.csv("UNgene_list.csv",header = T))
node[1]


node_used <- node
net_used <- net
k1 <- which(net_used[,1] %in% node_used)   
k2 <- which(net_used[,2] %in% node_used)   
length(intersect(k1,k2))   
used <- net_used[intersect(k1,k2),]


library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP, remove.loops = T, remove.multiple = T)  
ed <- as_edgelist(p1, names = TRUE)


clu <- components(p1)
groups(clu)
g <- p1

plot(g, layout=layout.fruchterman.reingold,  
     vertex.size=4,   
     vertex.label = V(g)$name,  
     vertex.label.cex=0.7,  
     vertex.label.dist=0.4,  
     vertex.label.color = "black"   
     )


component <- groups(clu)$'1'
# write.csv(component, file = "UNgene_component.csv", row.names = F)


node_used <- component
net_used <- net
k1 <- which(net_used[,1] %in% node_used)   
k2 <- which(net_used[,2] %in% node_used)   
length(intersect(k1,k2))  
used <- net_used[intersect(k1,k2),]


library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP, remove.loops = T, remove.multiple = T)  # 最终的数对
ed <- as_edgelist(p1, names = TRUE)
# write.csv(ed,"UNgene_comp_net.csv",row.names = F,quote = F)


g <- p1
plot(g, layout=layout.fruchterman.reingold, 
     vertex.size=4,  
     vertex.label = V(g)$name, 
     vertex.label.cex=0.7, 
     vertex.label.dist=0.4, 
     vertex.label.color = "black"  
)



# 20210913  UNgene_component - expression data ----------------------------------------------------------------

Data <- read.table("TCGA_pro_outcome_TN_log_UN901.txt",header=T,sep='\t', check.names = F)
dim(Data)   # 902   224
gene <- as.vector(rownames(Data)[-1])    


DE_outcome <- read.csv("UNgene_component.csv", header = T, sep=',')   
Degene <- DE_outcome[,1] 
dim(as.matrix(Degene))


data <- Data[Degene,]
all_data <- rbind(Data[1,], data)
dim(all_data)   


# scale -------------------------------------------------------------------

all_data_scale <- rbind(all_data[1,], t(scale(t(all_data[-1,]))))
dim(all_data_scale)    # 793 224
View(all_data_scale[,1:10])
# write.table(all_data_scale, file = "TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt",quote = F, sep = "\t")

