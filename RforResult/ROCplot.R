## clear
rm(list = ls())

## package
library(pROC)
library(ggplot2)
library(ROCR)  
## data
setwd('D:\\E\\...\\Data\\RTCGA\\result\\A_test')

roc1 <- read.csv("A_test_10780.csv", header=TRUE, sep = ',')
roc2 <- read.csv("A_test_42568.csv", header=TRUE, sep = ',')
roc3 <- read.csv("A_test_38959.csv", header=TRUE, sep = ',')
roc4 <- read.csv("A_test_20437.csv", header=TRUE, sep = ',')
roc5 <- read.csv("A_test_10797.csv", header=TRUE, sep = ',')
roc6 <- read.csv("A_test_21422.csv", header=TRUE, sep = ',')
roc7 <- read.csv("A_test_26910.csv", header=TRUE, sep = ',')
roc8 <- read.csv("A_test_45827.csv", header=TRUE, sep = ',')
roc9 <- read.csv("A_test_65194.csv", header=TRUE, sep = ',')


Roc1 <- roc(roc1[,2],roc1[,1])
g <- ggroc(Roc1)
g
g + theme_minimal() +  
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="red", 
               linetype=6) 
gl <- ggroc(Roc1, legacy.axes = TRUE)
gl

Roc2 <- roc(roc2[,2],roc2[,1])
Roc3 <- roc(roc3[,2],roc3[,1])
Roc4 <- roc(roc4[,2],roc4[,1])
Roc5 <- roc(roc5[,2],roc5[,1])
Roc6 <- roc(roc6[,2],roc6[,1])
Roc7 <- roc(roc7[,2],roc7[,1])
Roc8 <- roc(roc8[,2],roc8[,1])
Roc9 <- roc(roc9[,2],roc9[,1])

g2 <- ggroc(list(GSE10780=Roc1, 
                 GSE42568=Roc2, 
                 GSE38959=Roc3, 
                 GSE20437=Roc4,
                 GSE10797=Roc5,
                 GSE21422=Roc6, 
                 GSE26910=Roc7, 
                 GSE45827=Roc8,
                 GSE65194=Roc9),
            legacy.axes = TRUE)
g2


labels=c("GSE10780", "GSE42568", "GSE38959", "GSE20437", "GSE10797", "GSE21422", "GSE26910", "GSE45827", "GSE65194")

g2 + annotate(geom = "segment", 
              x = 0, y = 0, xend =1, yend = 1, 
              colour = "gray", size = 0.5) +
  scale_fill_discrete(labels) +
  theme_gray() + coord_equal() +
  theme(legend.position = c(0.70,0.35), #legend.position = 'inside',
        legend.text = element_text(color = 'black',size = 10),
        axis.text = element_text(color = 'black',size = 15),
        axis.text.x = element_text(angle = 0),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'))

