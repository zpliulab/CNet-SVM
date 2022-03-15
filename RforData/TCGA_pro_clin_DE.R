
library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)


# …Ë÷√¬∑æ∂ --------------------------------------------------------------------

setwd("D:\\E\\...Data\\RTCGA")


# BRCA and Normal -------------------------------------------------------

brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", FALSE) 
brca <- TCGAutils::splitAssays(brca, c('01','11'))
xdata.raw <- t(cbind(assay(brca[[1]]), assay(brca[[2]])))
dim(xdata.raw)    # 1205 20501



# Get matches between survival and assay data
class.v        <- TCGAbiospec(rownames(xdata.raw))$sample_definition %>% factor
names(class.v) <- rownames(xdata.raw)


xdata_raw <- xdata.raw %>%
{ (apply(., 2, sd) != 0) } %>%
{ xdata.raw[, .] }

dim(xdata_raw)    # 1205 20222

small_subset <- colnames(xdata.raw)
dim(small_subset)    # 1205 20222


xdata <- xdata_raw[, small_subset[small_subset %in% colnames(xdata_raw)]]
dim(xdata)    # 1205 20222
xdatat <- t(xdata)

ydata <- class.v
class(ydata)




## Normal
sample_N <- rownames(xdata)[which(ydata == "Solid Tissue Normal")]
sample_N1 <- sample_N %>% 
  as_tibble() %>% 
  mutate(sample_N = substr(sample_N, 1, 12)) 
## All
sample <- rownames(xdata)
sample1 <- sample %>% 
  as_tibble() %>% 
  mutate(sample = substr(sample, 1, 12)) 

## Tumor
lab <- which(as.matrix(sample1[,2]) %in% as.matrix(sample_N1[,2]))
xdata_TN <- xdata[lab,]
dim(xdata_TN)    # 224 20222

ydata_TN <- rbind(as.matrix(rep(c("Tumor"), 112)), as.matrix(rep(c("Normal"), 112)))
colnames(ydata_TN) <- c("outcome")


data_TN <- cbind(ydata_TN, xdata_TN)
dim(data_TN)    # 224 20223
# View(data_TN[,1:10])
# write.table(t(data_TN), "TCGA_pro_outcome_TN.txt",quote=F,sep="\t")



setwd("D:\\E\\...Data\\RTCGA")

library(DESeq2)
library(limma)
library(pasilla)



# DATA --------------------------------------------------------------------

data <- read.table("TCGA_pro_outcome_TN.txt",header=T,sep='\t', check.names = F)
dim(data)   # 20223   224


xdatat <- t(as.matrix(apply(as.matrix(data[-1,]),1,function(x) as.numeric(x))))
colnames(xdatat) <- colnames(data)
xdatat[1,2]
xdata[2,1]
ydata <- data[1,]


class(xdatat)
class(xdata)


# DESeq2  -------------------------------------------------

exprSet <- round(xdatat)
dim(exprSet)    # 20222  1205    20222   224


ydata_TN <- rbind(as.matrix(rep(c("Tumor"), 112)), as.matrix(rep(c("Normal"), 112)))
colnames(ydata_TN) <- c("outcome")

group_list <- as.factor(ydata_TN)
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
colnames(colData) <- c("outcome")
dim(colData)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ outcome)

dds2 <- DESeq(dds)  
resultsNames(dds2)


res <-  results(dds2, contrast=c("outcome","Tumor","Normal"))
summary(res) 
res_order <- res[order(res$padj),]
res_order <- as.data.frame(res_order)
# write.csv(res_order,file= "DEG_res_order_TN.csv")


## FDR
res1 <-  results(dds2, alpha = 0.01) 

diff_gene_deseq2 <- subset(res1, padj < 0.01 & abs(log2FoldChange) > 3.321928) #
# log2(9)=3.169925,10-3.321928
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
dim(diff_gene_deseq2)
# write.csv(diff_gene_deseq2,file= "DEG_Tumor_vs_Normal_10.csv")



vst_dat <- vst(dds2, blind = TRUE)
dat111 <- assay(vst_dat)
dim(dat111)    # 20222   224
View(dat111[,1:10])
# write.csv(dat111,file= "deseq_nor.csv")



data0 <- rbind(t(ydata_TN), dat111)
dim(data0)    # 20223   224
View(data0[,1:10])
data0[3,3]


ydata_TN <- rbind(as.matrix(rep(c("0"), 112)), as.matrix(rep(c("1"), 112)))
colnames(ydata_TN) <- c("outcome")
data0 <- rbind(t(ydata_TN), dat111)
dim(data0)    # 20223   224
View(data0[,1:10])
# write.table(data0,"TCGA_pro_outcome_TN_log.txt",quote=F,sep="\t") 



# DEGs symbol ---------------------------------------------------------

Data <- read.table("TCGA_pro_outcome_TN_log.txt",header=T,sep='\t', check.names = F)
dim(Data)   # 20223   224
# View(Data[,1:10])
gene <- as.vector(rownames(Data)[-1])    # 20222
View(gene)    


DE_outcome <- read.csv("DEG_Tumor_vs_Normal_10.csv", header = T, sep=',')    # 489
Degene <- DE_outcome[,1] # 489
# write.csv(Degene, file = "DEgene_list.csv", row.names = F)


intersect(gene, Degene)
data <- Data[Degene,]
dim(data)    # 1633 1080
View(data[,1:10])


all_data <- rbind(Data[1,], data)
dim(all_data)    # 2568 1080
# write.table(all_data, file = "TCGA_pro_outcome_TN_log_DE489.txt",quote = F, sep = "\t")


