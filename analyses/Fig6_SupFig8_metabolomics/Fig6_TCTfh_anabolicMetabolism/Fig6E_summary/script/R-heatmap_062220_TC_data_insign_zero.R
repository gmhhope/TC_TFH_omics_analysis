############# This script is intended to perform regular heatmap function ###########
library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(RColorBrewer)

############# Set up folder location and read tables ###################
HT <- read.csv("./input/TC_data/TC_invivo_sum_table.csv", header = T, row.names = 1,check.names = F) #
out_pdf_file = paste("./output/Heatmap_","ranged_2_","clustered_091420_padj.05insign_zero_disp_num",".pdf",sep = "") #

colnames(HT)

colnames(HT)[7:10]
HT_filt <- HT[HT$ttest_pval_B6_TfhvsB6_Tn < 0.05 | HT$ttest_pval_TC_TfhvsB6_Tfh < 0.05 |
                HT$ttest_pval_TC_TfhvsTC_Tn < 0.05 | HT$ttest_pval_TC_TnvsB6_Tn < 0.05, ]

HT_filt[HT_filt$ttest_pval_B6_TfhvsB6_Tn > 0.05, "log2FC_B6_TfhvsTn"] = 0
HT_filt[HT_filt$ttest_pval_TC_TfhvsB6_Tfh > 0.05, "log2FC_Tfh_TCvsB6"] = 0
HT_filt[HT_filt$ttest_pval_TC_TfhvsTC_Tn > 0.05, "log2FC_TC_TfhvsTn"] = 0
HT_filt[HT_filt$ttest_pval_TC_TnvsB6_Tn > 0.05, "log2FC_Tn_TCvsB6"] = 0


data = HT_filt[,7:10]

data[data > 2] =2
data[data < -2] = -2
############# Z-score tranformation (optional)  #################

# zscore <- function(x) {
#   z <- (x - rowMeans(x))/ apply(x,1,sd)
#   return(z)
# }
# 
# z_transformed <-zscore(data)


annotation_col = data.frame(Var1 = colnames(HT)[7:10]) #c(rep("Tact",4),rep("Tfh",3),rep("Th1",4),rep("Th17",4)
rownames(annotation_col) = colnames(data)

annotation_row = row.names(HT_filt) 
###############   Write your output file name   ###################
pdf(file=out_pdf_file, width =20,height=40, paper = "special",onefile=FALSE)


### Either z_transformed or HT ###
paletteLength <- 50
myColor <- colorRampPalette(c("yellow", "white", "purple"))(paletteLength)

myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1),  #min(data)
              seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))


pheatmap(data, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col, clustering_method = "complete", cellwidth = 19.5,cellheight = 10,
         color = myColor, breaks = myBreaks, display_numbers = TRUE)
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()
