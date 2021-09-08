############# This script is intended to perform regular heatmap function ###########
library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)

############# Set up folder location and read tables ###################
setwd("/Volumes/LM_MG_drive_one/2016Immunity_Metabolomics/X6_cluster_analysis/")
HT <- read.csv("combined_df_only_TfhvsTact_mod.csv", header = T, row.names = 2,check.names = F)
colnames(HT)

colnames(HT)[5:11]
df_l = list()

HT_filt <- HT[HT$padj_lim_ttest_tscore__TactvsTfh < 0.05, ]
HT_filt <- HT_filt[order(HT_filt$pathway,HT_filt$ttest_tscore__TactvsTfh),]
data = HT_filt[,5:11]
############# Z-score tranformation (optional)  #################

zscore <- function(x) {
  z <- (x - rowMeans(x))/ apply(x,1,sd)
  return(z)
}

z_transformed <- t(scale(t(data)))


annotation_col = data.frame(Var1 = c(rep("Tact",4),rep("Tfh",3)))
rownames(annotation_col) = colnames(data)

annotation_row = data.frame(Var2 = HT_filt$pathway) 
rownames(annotation_row) = rownames(HT_filt)

###############   Write your output file name   ###################
out_pdf_file = paste("Heatmap_","ordered_only_TfhvsTact_trial","clustered_06122020_padj.05",".pdf",sep = "")
pdf(file=out_pdf_file, width =20,height=40, paper = "special",onefile=FALSE)


### Either z_transformed or HT ###
pheatmap(z_transformed, cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = annotation_col, annotation_row = annotation_row, clustering_method = "complete", cellwidth = 10,cellheight = 10,
         color = colorRampPalette(c("blue", "white", "red"))(50)) 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()
