
#  ˝æ› ‰»Î --------------------------------------------------------------------

setwd("D:\\E\\...Data\\RTCGA")

Data = read.table("TCGA_pro_outcome_TN_log.txt", header = T, check.names = FALSE)
dim(Data)    # 16145   224
gene <- as.vector(rownames(Data)[-c(1)])


DEData = read.table("TCGA_pro_outcome_TN_log_DE489.txt", header = T, check.names = FALSE)
dim(DEData)    # 391 1080
# View(DEData[,1:10])
Degene <-  as.vector(rownames(DEData)[-c(1)])

intersect(gene, Degene)
sub <- setdiff(gene, Degene)     


# biomarker ---------------------------------------------------------------

setwd('D:\\E\\...Data')


malacards <- as.matrix(read.csv("malacards82.csv", header = T, sep=','))
GEDFN <- as.matrix(read.csv("GEDFN169.csv", header = T, sep=','))
mama_70 <- as.matrix(read.csv("mamaprint70.csv", header = T, sep=',') )
KEGG_147 <- as.matrix(read.csv("KEGG147.csv", header = T, sep=',') )
tf <- as.matrix(read.csv("tfgene119.csv", header = T, sep = ','))


# union  ------------------------------------------------------------

add_gene <- as.matrix( union( union(union(union(intersect(sub, malacards), 
                                                intersect(sub, GEDFN)), 
                                          intersect(sub, mama_70)), 
                                    intersect(sub, KEGG_147)), 
                              intersect(sub, tf) ) 
                       )

rownames(add_gene) <- add_gene[,1] 

add_data <- Data[rownames(add_gene),]
dim(add_data)    # 412  224

all_data <- rbind(DEData, add_data)
dim(all_data)    # 902 224
View(all_data[,1:10])

geneall <- rownames(all_data)[-1]
View(geneall)

# write.csv(geneall, file = 'UNgene_list.csv', row.names = F)
# write.table(all_data, file = "TCGA_pro_outcome_TN_log_UN901.txt",quote = F, sep = "\t")


