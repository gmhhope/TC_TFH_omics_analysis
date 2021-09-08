############# This script is intended to perform regular heatmap function ###########
library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(rlist)
library("org.Mm.eg.db") #BiocManager::install("org.Mm.eg.db")

############# Set up folder location and read tables ###################
HT <- read.csv("input/full_RNAseq_LRT_added_mean_added_annot_meta_stat_v3.csv", header = T, check.names = F, stringsAsFactors = FALSE)

treatment_grp = c("Tfh_B6","Tn_B6","Tfh_TC","Tn_TC")
num_grp = c(4,4,3,3)
sample_num = sum(num_grp)
padj_thre = 0.05

colnames(HT)[1] <- "X"
HT = HT[!duplicated(HT$X),]
row.names(HT) <- HT$X
HT <- HT[HT$padj_LRT < padj_thre & !is.na(HT$padj_LRT),] #

#genelist.df = read.csv("input/Tfh_specific_genes_Immunity2014_v2.csv", header = T, check.names = F)
colnames(HT)[6:19]# %in% genelist.df$X

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


data = HT[,6:19]
write.csv(data,"output/data.csv")


############# scale  #################

scaled <-scale(t(data),center = TRUE,scale = TRUE)
scaled <- t(scaled)
write.csv(scaled,"output/scaled.csv")

annotation_col = data.frame(Var1 = c(rep(treatment_grp[1],num_grp[1]),
                                     rep(treatment_grp[2],num_grp[2]),
                                     rep(treatment_grp[3],num_grp[3]),
                                     rep(treatment_grp[4],num_grp[4])))
rownames(annotation_col) = colnames(data)
colnames(HT)
annotation_row = HT[,c('annot','stat_annot','stat_interaction','metabolism')]
###############   Write your output file name   ###################
out_pdf_file = paste("output/Heatmap_","all_saoutmples_","quantiled_09112020_selected_padj05",".pdf",sep = "")
pdf(file=out_pdf_file, width =6,height = 6, paper = "special",onefile=FALSE)



range(scaled) #check the range of z_transformed
breaksList = seq(min(scaled), max(scaled), length.out = 10) #seq(min(mat), max(mat), length.out = 10)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(scaled, n = 11)
cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(mat_breaks)) #rev()


### Either z_transformed or HT ###
pheatmap(scaled, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = TRUE, 
         color = cell_colors,
         breaks = mat_breaks,
         annotation_col = annotation_col, 
         annotation_row = annotation_row, 
         clustering_method = "complete") 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()