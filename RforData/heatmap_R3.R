## Mar. 2023, Followed "TCGA_pro_clin_DE_R3.R"
## heatmap_R3.R, is used to see the selected 32 biomarkers's performance on whole 1205 sample.


## clear
rm(list = ls())

filepath <- '/Users/lilingyu/E/PhD/'
# filepath <- '/home/lly/'

# load data ---------------------------------------------------------------

setwd(paste(filepath,'R/SVM/Data/RTCGA',sep = ""))
x0 = read.table("TCGA_pro_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names=F)
dim(x0)    # 20223  1205
View(x0[1:10,1:10])



# Extract 32 biomarker gene symbol ---------------------------------------------------------

Data <- read.table("TCGA_pro_outcome_TN_log_R3.txt",header=T,sep='\t', check.names = F)
dim(Data)   # 20223  1205
gene <- as.vector(rownames(Data)[-1])    # 20222
View(gene)    


DE_outcome <- read.csv("result/featurenew/CsvdataFeatureSeleOnce/CNetSVM.csv", header = T, sep=',')
Degene <- DE_outcome[,1]
data <- Data[Degene,]


all_data <- rbind(Data[1,], data)
dim(all_data)
# write.table(all_data, file = "TCGA_pro_outcome_TN_log_DE489_R3_32biomarker.txt",quote = F, sep = "\t")


# Load data ---------------------------------------------------------------
x0 = read.table("TCGA_pro_outcome_TN_log_DE489_R3_32biomarker.txt", header = T, 
                sep='\t', fill=TRUE, strip.white = T, check.names=F)
rownames(x0)
x1 <- data.frame(t(x0))
colnames(x1) <- rownames(x0)
x2 <- x1[,-1]

View(rownames(x2))
View(x1[,1])
label <- cbind(rownames(x2),x1[,1])
colnames(label) <- c("sample","class")
# write.table(label, file = "phe_zhSPTB_R3.txt",quote=F,sep="\t",row.names = F)


# heatmap -----------------------------------------------------------------
library(pheatmap)

selected <- t(x2)
Label = read.table("phe_zhSPTB_R3.txt",header = TRUE, sep = "\t")
Label <- factor(Label[,"class"])
Label <- data.frame(Label)
rownames(Label) = colnames(selected)


selectednum <- apply(selected, 2, as.numeric)
matrix.smooth <- apply(selectednum,1,function(x){smooth.spline(x)$y})
normalized_matrix <- apply(matrix.smooth,1,function(x){(x-mean(x))/sd(x)})
pheatmap (normalized_matrix, show_colnames=F,cluster_row=F, cluster_col=F)


colnames(normalized_matrix) <- colnames(selectednum) <- colnames(selected)
rownames(normalized_matrix) <- rownames(selectednum) <- rownames(selected)
pp <- pheatmap(normalized_matrix,    # normalized_matrix    selectednum
               annotation_col = Label, 
               # color = colorRampPalette(c("blue", "white","red"))(100),
               fontsize_row = 8,
               # scale = "row", 
               # cutree_cols = 2,
               cluster_col = F,
               border_color = NA,
               cluster_row = F,
               show_colnames = F,
               show_rownames = T)

pp

