############# This script is intended to perform regular heatmap function ###########
library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(RColorBrewer)

############# Set up folder location and read tables ###################
HT <- read.csv("./input/TC_data/TC_invivo_sum_table.csv", header = T, check.names = F) #
metab_list = read.csv("./input/TC_data/metabolite_list_presented_from_TCinvivo.csv", stringsAsFactors = FALSE)
out_pdf_file = paste("./output/Heatmap_","ranged_2_","clustered_091420_padj.05insign_zero_disp_num_selected_ordered",".pdf",sep = "") #
metab_list$metabolite <- factor(metab_list$metabolite, levels = metab_list$metabolite)
HT <- merge(HT,metab_list, by.x = "BinBase name", by.y = "metabolite", how = "inner")
row.names(HT) <- HT$`BinBase name`
HT$`BinBase name` <- factor(HT$`BinBase name`, levels = metab_list$metabolite)
HT <- HT[order(HT$`BinBase name`),]
colnames(HT)

colnames(HT)[8:11]
HT_filt <- HT[HT$ttest_pval_B6_TfhvsB6_Tn < 0.05 | HT$ttest_pval_TC_TfhvsB6_Tfh < 0.05 |
                HT$ttest_pval_TC_TfhvsTC_Tn < 0.05 | HT$ttest_pval_TC_TnvsB6_Tn < 0.05, ]

HT_filt[HT_filt$ttest_pval_B6_TfhvsB6_Tn > 0.05, "log2FC_B6_TfhvsTn"] = 0
HT_filt[HT_filt$ttest_pval_TC_TfhvsB6_Tfh > 0.05, "log2FC_Tfh_TCvsB6"] = 0
HT_filt[HT_filt$ttest_pval_TC_TfhvsTC_Tn > 0.05, "log2FC_TC_TfhvsTn"] = 0
HT_filt[HT_filt$ttest_pval_TC_TnvsB6_Tn > 0.05, "log2FC_Tn_TCvsB6"] = 0


data = HT_filt[,8:11]


data[data > 2] =2
data[data < -2] = -2
############# Z-score tranformation (optional)  #################

# zscore <- function(x) {
#   z <- (x - rowMeans(x))/ apply(x,1,sd)
#   return(z)
# }
# 
# z_transformed <-zscore(data)


annotation_col = data.frame(Var1 = colnames(HT)[8:11]) #c(rep("Tact",4),rep("Tfh",3),rep("Th1",4),rep("Th17",4)
rownames(annotation_col) = colnames(data)

annotation_row = row.names(HT_filt) 
###############   Write your output file name   ###################
pdf(file=out_pdf_file, width =20,height=40, paper = "special",onefile=FALSE)


### Either z_transformed or HT ###
paletteLength <- 50
myColor <- colorRampPalette(c("yellow", "white", "purple"))(paletteLength)

myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1),  #min(data)
              seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))


pheatmap(data, cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE,
         annotation_col = annotation_col, clustering_method = "complete", cellwidth = 19.5,cellheight = 10,
         color = myColor, breaks = myBreaks, display_numbers = TRUE)
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()
