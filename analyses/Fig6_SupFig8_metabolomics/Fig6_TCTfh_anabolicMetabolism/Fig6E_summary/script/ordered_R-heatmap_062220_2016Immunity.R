############# This script is intended to perform regular heatmap function ###########
library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(RColorBrewer)

############# Set up folder location and read tables ###################
HT <- read.csv("./input/2016Immunity/2016Immunity_TfhvsTact_sum_table.csv", header = T, check.names = F)
out_pdf_file = paste("./output/Heatmap_","ranged_2_","2016_Immunity_clustered_091420_padj.05_insign_zero_display_num_selected",".pdf",sep = "")
metab_list = read.csv("./input/2016Immunity/metabolite_list_presented_from_2016Immunity.csv")
HT <- merge(HT,metab_list, by = "metabolite", how = "inner")
HT$metabolite = factor(HT$metabolite, levels = metab_list$metabolite)
HT <- HT[order(HT$metabolite),]
colnames(HT)[6]
HT_filt <- HT[HT$ttest_pval_TactvsTfh < 0.05, ]
HT[HT$ttest_pval_TactvsTfh < 0.05, "log2FC_TfhvsTact"] = 0
data = as.data.frame(HT_filt[,6])
colnames(data) = colnames(HT)[6]
row.names(data) = HT$metabolite

data[data > 2] = 2
data[data < -2] = -2
############# Z-score tranformation (optional)  #################

# zscore <- function(x) {
#   z <- (x - rowMeans(x))/ apply(x,1,sd)
#   return(z)
# }
# 
# z_transformed <-zscore(data)


annotation_col = data.frame(Var1 = colnames(HT)[6]) #c(rep("Tact",4),rep("Tfh",3),rep("Th1",4),rep("Th17",4)
rownames(annotation_col) = colnames(data)

annotation_row = row.names(HT_filt) 
###############   Write your output file name   ###################
pdf(file=out_pdf_file, width =20,height=40, paper = "special",onefile=FALSE)


### Either z_transformed or HT ###
paletteLength <- 50
myColor <- colorRampPalette(c("yellow", "white", "purple"))(paletteLength)

myBreaks <- c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))


pheatmap(data, cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE,
         annotation_col = annotation_col, clustering_method = "complete", cellwidth = 19.5,cellheight = 10,
         color = myColor, breaks = myBreaks, display_numbers = TRUE)
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()
