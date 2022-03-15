library(dplyr)        
library(tidyr)
library(tidyverse)   



# Find gene expression data based on feature --------------------------------------------------
rm(list = ls())

# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE76275')
# Data1 = read.table("GSE76275_outcome_scale.txt", header = T, check.names = FALSE)


# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE15852')
# Data1 = read.table("GSE15852_outcome_scale.txt", header = T, check.names = FALSE)


# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE25407')
# Data1 = read.table("GSE25407_outcome_scale.txt", header = T, check.names = FALSE)


# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE10180')
# Data1 = read.table("GSE10180_outcome_scale.txt", header = T, check.names = FALSE)


# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE10780')
# Data1 = read.table("GSE10780_outcome_scale.txt", header = T, check.names = FALSE)


# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE42568')
# Data1 = read.table("GSE42568_outcome_scale.txt", header = T, check.names = FALSE)


# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE70905')
# Data1 = read.table("GSE70905_outcome_scale.txt", header = T, check.names = FALSE)


# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE10797')
# Data1 = read.table("GSE10797_outcome_scale.txt", header = T, check.names = FALSE)


# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE61304')
# Data1 = read.table("GSE61304_outcome_scale.txt", header = T, check.names = FALSE)

# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE45827')
# Data1 = read.table("GSE45827_outcome_scale.txt", header = T, check.names = FALSE)

# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE65194')
# Data1 = read.table("GSE65194_outcome_scale.txt", header = T, check.names = FALSE)

# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE38959')
# Data1 = read.table("GSE38959_outcome_scale.txt", header = T, check.names = FALSE)

# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE36693')
# Data1 = read.table("GSE36693_outcome_scale.txt", header = T, check.names = FALSE)

# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE21422')
# Data1 = read.table("GSE21422_outcome_scale.txt", header = T, check.names = FALSE)

# setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE26910')
# Data1 = read.table("GSE26910_outcome_scale.txt", header = T, check.names = FALSE)

setwd('D:\\E\\博士\\R_程序\\SVM\\Data\\GEO\\GSE20437')
Data1 = read.table("GSE20437_outcome_scale.txt", header = T, check.names = FALSE)


dim(Data1)    # 21836   265

setwd('D:\\...\\Data\\RTCGA\\result')


myfile = list.files("feature")                 
dir = paste("./feature/", myfile, sep = "")     
n = length(dir)                                  
for (i in 1:n) {
  gene = read.csv(file = dir[i],
                  header = T,
                  sep = ",")
  
  
  colnames(gene) <- c('gene')
  Data2 <- cbind(rownames(Data1), Data1)
  colnames(Data2) <- c('gene', colnames(Data1))
  
  genedata <- merge(gene, Data2, by = "gene")
  genedata1 <-
    genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(Data1[1,], genedata1)
  rownames(genedata2) <- c('Lable', rownames(genedata1))
  
  # i <- 4
  name <- dir[i]
  name1 <- str_split_fixed(name, "./", 2)
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[.]", 2)
  name4 <- name3[1]
  # name5 <- str_c("Feature76275_", name4)
  # name5 <- str_c("Feature15852_", name4)
  # name5 <- str_c("Feature25407_", name4)
  # name5 <- str_c("Feature10810_", name4)
  # name5 <- str_c("Feature10780_", name4)
  # name5 <- str_c("Feature42568_", name4)
  # name5 <- str_c("Feature70905_", name4)
  # name5 <- str_c("Feature10797_", name4)
  # name5 <- str_c("Feature61304_", name4)
  # name5 <- str_c("Feature45827_", name4)
  # name5 <- str_c("Feature65194_", name4)
  # name5 <- str_c("Feature38959_", name4)
  # name5 <- str_c("Feature36693_", name4)
  # name5 <- str_c("Feature21422_", name4)
  # name5 <- str_c("Feature26910_", name4)
  name5 <- str_c("Feature20437_", name4)
  
  
  path <-
    paste("./featuredatanew/", paste(name5, ".txt", sep = ""), sep = "")
  write.table(genedata2, path, quote = F, sep = "\t")
  
  
  # setwd('D:\\...\\Data\\RTCGA')
  data1 = read.table(
    "D:\\E\\...\\Data\\RTCGA\\TCGA_pro_outcome_TN_log_comp_UNgene_scale.txt",
    header = T,
    check.names = FALSE
  )
  dim(data1)    # 793 224
  gene <- as.matrix(genedata[, 1])
  
  colnames(gene) <- c('gene')
  data2 <- cbind(rownames(data1), data1)
  colnames(data2) <- c('gene', colnames(data1))

  
  genedata <- merge(gene, data2, by = "gene")
  genedata1 <-
    genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(data1[1,], genedata1)
  rownames(genedata2) <- c('Lable', rownames(genedata1))
  
  
  # name5 <- str_c("Feature_tcga_76275_", name4)
  # name5 <- str_c("Feature_tcga_15852_", name4)
  # name5 <- str_c("Feature_tcga_25407_", name4)
  # name5 <- str_c("Feature_tcga_10810_", name4)
  # name5 <- str_c("Feature_tcga_10780_", name4)
  # name5 <- str_c("Feature_tcga_42568_", name4)
  # name5 <- str_c("Feature_tcga_70905_",name4)
  # name5 <- str_c("Feature_tcga_10797_", name4)
  # name5 <- str_c("Feature_tcga_61304_", name4)
  # name5 <- str_c("Feature_tcga_45827_", name4)
  # name5 <- str_c("Feature_tcga_65194_", name4)
  # name5 <- str_c("Feature_tcga_38959_", name4)
  # name5 <- str_c("Feature_tcga_36693_", name4)
  # name5 <- str_c("Feature_tcga_21422_", name4)
  # name5 <- str_c("Feature_tcga_26910_", name4)
  name5 <- str_c("Feature_tcga_20437_", name4)
  
  
  path <- paste("./featuredatanew/", paste(name5,".txt", sep=""), sep="")
  write.table(genedata2, path, quote = F, sep="\t")
  
}              

