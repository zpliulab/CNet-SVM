## Mar. 2023, New version of Rcode "TCGA_pro_clin_DE.R"
## TCGA_pro_clin_DE_R3.R, it accesses the whole data of breast cancer, 1205 samples.
## TCGA_pro_clin_DE.R", it only accesses 112 normal and 112 tumor samples.


# R package --------------------------------------------------------------------

# if (!require("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("curatedTCGAData")


# R package --------------------------------------------------------------------

library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)


# load data --------------------------------------------------------------------

setwd("D:\\E\\博士\\R_程序\\SVM\\Data\\RTCGA")


# Tumor and Normal samples -------------------------------------------------------

brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", FALSE)
brca <- TCGAutils::splitAssays(brca, c('01','11'))
xdata.raw <- t(cbind(assay(brca[[1]]), assay(brca[[2]])))

dim(xdata.raw)    # 1205 20501
# View(xdata.raw[,1:10])


# Get matches between survival and assay data
class.v        <- TCGAbiospec(rownames(xdata.raw))$sample_definition %>% factor
names(class.v) <- rownames(xdata.raw)


xdata_raw <- xdata.raw %>%
{ (apply(., 2, sd) != 0) } %>%
{ xdata.raw[, .] }


dim(xdata_raw)    # 1205 20222
# View(xdata_raw[,1:10])


small_subset <- colnames(xdata.raw)
dim(small_subset)    # 1205 20222
# View(small_subset[,1:10])

xdata <- xdata_raw[, small_subset[small_subset %in% colnames(xdata_raw)]]
dim(xdata)    # 1205 20222
# View(xdata[,1:10])
xdatat <- t(xdata)

ydata <- class.v
# View(ydata)
class(ydata)



# Use the whole dataset -------------------------------------------------------------

xydata <- cbind(ydata, xdata)
# dim(xydata)    # 1205 20223
# View(xydata[,1:10])
xydata[which(xydata[,1] == "1"), 1] <- c("Tumor")
xydata[which(xydata[,1] == "2"), 1] <- c("Normal")
# View(xydata[,1:10])

# write.table(t(xydata), "TCGA_pro_outcome.txt",quote=F,sep="\t")


setwd("D:\\E\\博士\\R_程序\\SVM\\Data\\RTCGA")

# DEGs packages --------------------------------------------------------------------


# BiocManager::install("DESeq2")
library(DESeq2)
library(limma)
library(pasilla)



# Load Data --------------------------------------------------------------------
data <- read.table("TCGA_pro_outcome.txt",header=T,sep='\t', check.names = F)
dim(data)   # 20223  1205
# View(data[,1:10])
# allgene <- rownames(data)[-1]
# View(allgene)
# write.csv(allgene, file = "allgene_list.csv", row.names = F)

xdatat <- t(as.matrix(apply(as.matrix(data[-1,]),1,function(x) as.numeric(x))))
colnames(xdatat) <- colnames(data)
# View(xdatat[,1:10])
xdatat[1,2]
xdata[2,1]
ydata <- data[1,]


# DESeq2 performe DEG analysis -------------------------------------------------

exprSet <- round(xdatat)
dim(exprSet)    # 20222  1205
# View(exprSet[,1:10])


ydata_TN <- as.matrix(as.character(t(ydata)))
colnames(ydata_TN) <- c("outcome")
# View(ydata_TN)

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
# plotMA(res)
# 
# # BiocManager::install("apeglm")
# library(apeglm)
# resLFC <- lfcShrink(dds2, coef="outcome_Tumor_vs_Normal", type="apeglm")
# plotMA(resLFC) 

res_order <- res[order(res$padj),]
res_order <- as.data.frame(res_order)
# write.csv(res_order,file= "DEG_res_order_TN_R3.csv")



# Set threshold -----------------------------------------------------------

## FDR adjust
res1 <-  results(dds2, alpha = 0.01)
# write.csv(res1,file= "DEG_res_R3.csv")


diff_gene_deseq2 <- subset(res1, padj < 0.01 & abs(log2FoldChange) > 3.321928) #
# log2(9)=3.169925,10-3.321928
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
dim(diff_gene_deseq2)
# write.csv(diff_gene_deseq2,file= "DEG_Tumor_vs_Normal_10_R3.csv")



# Normalization -------------------------------------------------------------------

# normalized_counts <- counts(dds2, normalized=TRUE)
# View(normalized_counts[,1:10])
vst_dat <- vst(dds2, blind = TRUE)
dat111 <- assay(vst_dat)
dim(dat111)    # 20222   224
View(dat111[,1:10])
# write.csv(dat111,file= "deseq_nor_R3.csv")


data0 <- rbind(t(ydata_TN), dat111)
dim(data0)    # 20223   224
View(data0[,1:10])
data0[3,3]


colnames(ydata_TN) <- c("outcome")
data0 <- rbind(t(ydata_TN), dat111)
dim(data0)    # 20223   224
View(data0[,1:10])
# write.table(data0,"TCGA_pro_outcome_TN_log_R3.txt",quote=F,sep="\t") 
