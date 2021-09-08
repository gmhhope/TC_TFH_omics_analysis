############# This script is intended to perform regular heatmap function ###########
library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(rlist)
library("org.Mm.eg.db") #BiocManager::install("org.Mm.eg.db")

############# Set up folder location and read tables ###################
#setwd("")
HT <- read.csv("input/either_R_N_OR_M_N_cleaned_combn4further_analyses_padj0.1.csv", header = T, check.names = F)

treatment_grp = c("Tfh_Nanostring","Tfh_mArray","Tfh_RNAseq")
num_grp = c(1,1,1)
sample_num = sum(num_grp)
padj_thre = 0.05

colnames(HT)[2] <- "X"
HT = HT[!duplicated(HT$X),]
row.names(HT) <- HT$X
# HT <- HT[HT$RNAseq_padj < padj_thre & !is.na(HT$RNAseq_padj)
#          & HT$mArray_padj < padj_thre & !is.na(HT$RNAseq_padj),] #

#genelist.df = read.csv("input/Tfh_specific_genes_Immunity2014_v2.csv", header = T, check.names = F)
colnames(HT)[c(2,6,7,8,1)]# %in% genelist.df$X

# if(T) {
#   HT_filt$GeneName <- mapIds(org.Mm.eg.db,
#                              keys=sapply(rownames(HT_filt),function(x)as.character(x)),
#                              column="GENENAME",
#                              keytype="SYMBOL",
#                              multiVals="first")
#   for (i in 1:length(row.names(HT_filt))) {
#     row.names(HT_filt)[i] <- paste(row.names(HT_filt)[i], "_", HT_filt$GeneName[i], sep = "")
#   }
# }
write.csv(HT,"output/HT.csv")


data = HT[,c(6,7,8)]
write.csv(data,"output/data.csv")


############# scale  #################
scaled = data.matrix(data)
#scaled <-scale(t(data),center = TRUE,scale = TRUE)
#scaled <- t(scaled)
#write.csv(scaled,"output/scaled.csv")

annotation_col = data.frame(Var1 = c(rep(treatment_grp[1],num_grp[1]),
                                     rep(treatment_grp[2],num_grp[2]),
                                     rep(treatment_grp[3],num_grp[3])))
rownames(annotation_col) = colnames(data)
colnames(HT)
annotation_row = data.frame("pathway" = HT[,c("pathway")])
rownames(annotation_row) = HT[,"X"]

###############   Write your output file name   ###################
out_pdf_file = paste("output/Heatmap_","","quantiled_01302020_selected_padj0.1",".pdf",sep = "")
pdf(file=out_pdf_file, width =5,height = 10, paper = "special",onefile=FALSE)

#regular break
if(F) {
  range(scaled) #check the range of z_transformed
  breaksList = seq(-2, 2, 0.1)
  cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(breaksList)) #rev()
  mat_breaks <- seq(min(scaled), max(scaled), length.out = 0.1)
  
}

range(scaled) #check the range of z_transformed
breaksList = seq(min(scaled), max(scaled), length.out = 10) #seq(min(mat), max(mat), length.out = 10)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n)) 
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(scaled, n = 11) #The quantile function require matrix rather than data.frame! You may get error!
cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(mat_breaks)) #rev()


### Either z_transformed or HT ###
pheatmap(scaled, 
         display_numbers = TRUE,
         cluster_rows = FALSE, 
         show_rownames = TRUE, 
         cluster_cols = TRUE, 
         color = cell_colors,
         breaks = mat_breaks,
         annotation_col = annotation_col, 
         annotation_row = annotation_row, 
         clustering_method = "complete") 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()