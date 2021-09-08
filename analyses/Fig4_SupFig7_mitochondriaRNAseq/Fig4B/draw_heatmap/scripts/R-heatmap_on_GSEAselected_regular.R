############# This script is intended to perform regular heatmap function ###########
library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(rlist)
library("org.Mm.eg.db") #BiocManager::install("org.Mm.eg.db")

############# Set up folder location and read tables ###################
#setwd("")
HT <- read.csv("input/ROS_selected.csv", header = T, row.names = 1,check.names = F)
out_pdf_file = paste("output/Heatmap_","ROS_selected_","RNAseq",".pdf",sep = "")

treatment_grp = c("log2FC_Tn_TCvsB6","log2FC_Tfh_TCvsB6","log2FC_B6_TfhvsTn","log2FC_TC_TfhsTn")
num_grp = c(1,1,1,1)
sample_num = sum(num_grp)
#padj_thre = 0.05

# colnames(HT)[1] <- "X"
# HT = HT[!duplicated(HT$X),]
# row.names(HT) <- HT$X
# HT <- HT[HT$RNAseq_pdj < padj_thre & !is.na(HT$RNAseq_padj)
#          & HT$mArray_padj < padj_thre & !is.na(HT$RNAseq_padj),] #

#genelist.df = read.csv("input/Tfh_specific_genes_Immunity2014_v2.csv", header = T, check.names = F)
colnames(HT)[c(5:8)]# %in% genelist.df$X

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


data = HT[,c(5:8)]
write.csv(data,"output/data.csv")


############# scale  #################
scaled = data.matrix(data)
#scaled <-scale(t(data),center = TRUE,scale = TRUE)
#scaled <- t(scaled)
#write.csv(scaled,"output/scaled.csv")

annotation_col = data.frame(Var1 = c(rep(treatment_grp[1],num_grp[1]),
                                     rep(treatment_grp[2],num_grp[2]),
                                     rep(treatment_grp[3],num_grp[3]),
                                     rep(treatment_grp[4],num_grp[4])))
rownames(annotation_col) = colnames(data)
colnames(HT)
# annotation_row = data.frame("pathway" = HT[,c("pathway")])
# rownames(annotation_row) = HT[,"X"]

###############   Write your output file name   ###################
pdf(file=out_pdf_file,  paper = "special",onefile=FALSE, width =5,height = 30) #width =5,height = 10,

#regular break
if(T) {
  range(scaled) #check the range of z_transformed
  paletteLength <- 50
  cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(paletteLength) #rev()
  #cell_colors <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  #mat_breaks <- c(seq(min(scaled), 0, length.out=ceiling(paletteLength/2) + 1), 
  #                seq(max(scaled)/paletteLength, max(scaled), length.out=floor(paletteLength/2)))
  mat_breaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(2/paletteLength, 2, length.out=floor(paletteLength/2)))
  
}

if (F) {
  range(scaled) #check the range of z_transformed
  breaksList = seq(min(scaled), max(scaled), length.out = 10) #seq(min(mat), max(mat), length.out = 10)

  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
    mat_breaks <- quantile_breaks(scaled, n = 11) #The quantile function require matrix rather than data.frame! You may get error!
    
  }
  cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(mat_breaks)) #rev()
  
}



### Either z_transformed or HT ###
pheatmap(scaled, 
         cellwidth = 24,
         cellheight = 12,
         display_numbers = TRUE,
         cluster_rows = TRUE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE, 
         color = cell_colors,
         breaks = mat_breaks,
         annotation_col = annotation_col, 
         #annotation_row = annotation_row, 
         clustering_method = "complete") 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()