library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(rlist)
library(gridExtra)
library(ggplot2)
library("org.Mm.eg.db") #BiocManager::install("org.Mm.eg.db")

HT <- read.csv("input/full_RNAseq_LRT_added_samples_removed.csv", header = T,check.names = F) # Input of the file, which basically remove all the columns with regularized transformed values
gene_df <- read.csv("input/FACS_RNAseq_gene_name_conv_v2.csv",header = T, stringsAsFactors = TRUE) # Handling the discrepancy of name convention between RNA-seq and FACS data
HT <- HT[HT$MouseSymbol %in% gene_df$RNAseq,]
show_trend = "not_show_trend" # not showing the genes that are not significant
row.names(HT) = HT$MouseSymbol

category_list <- unique(gene_df$Category) 

colnames(HT)[c(18,14,6,2)]# check whether the data columns in this area

colnames(HT)

if(show_trend == "not_show_trend") { # put all the non-significant values into 0, later n.s. are manually replace in the plot.
  HT[HT$padj_TCTfhvsTCTn > 0.05,"log2FoldChange_TCTfhvsTCTn"] = 0 
  HT[HT$padj_B6TfhvsB6Tn > 0.05,"log2FoldChange_B6TfhvsB6Tn"] = 0 
  HT[HT$padj_TCTfhvsB6Tfh > 0.05,"log2FoldChange_TCTfhvsB6Tfh"] = 0 
  HT[HT$padj_TCTnvsB6Tn > 0.05,"log2FoldChange_TCTnvsB6Tn"] = 0 
}

data = HT[,c(18,14,6,2)]


#annotation of the columns
treatment_grp = c("log2FC_Tn_TCvsB6","log2FC_Tfh_TCvsB6","log2FC_B6_TfhvsTn","log2FC_TC_TfhsTn")
num_grp = c(1,1,1,1)
sample_num = sum(num_grp)

annotation_col = data.frame(Var1 = c(rep(treatment_grp[1],num_grp[1]),
                                     rep(treatment_grp[2],num_grp[2]),
                                     rep(treatment_grp[3],num_grp[3]),
                                     rep(treatment_grp[4],num_grp[4])))
rownames(annotation_col) = colnames(data)

# plot multiple files #

plot_list=list()
for (a in category_list){
  #select data
  scaled = data.matrix(data[row.names(data) %in% gene_df[gene_df$Category==a,"RNAseq"],]) #
  scaled = scaled[complete.cases(scaled),]
  #regular break
  if(T) {
    range(scaled) #check the range of z_transformed
    paletteLength <- 50
    cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(paletteLength) #rev()
    #cell_colors <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    mat_breaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1), 
                    seq(2/paletteLength, 2, length.out=floor(paletteLength/2)))
    
  }
  
  x=pheatmap(scaled, 
             cellwidth = 24,
             cellheight = 12,
             display_numbers = TRUE,
             cluster_rows = TRUE, 
             show_rownames = TRUE, 
             show_colnames = FALSE,
             cluster_cols = FALSE, 
             color = cell_colors,
             breaks = mat_breaks,
             annotation_col = annotation_col, 
             main = a,
             #annotation_row = annotation_row, 
             clustering_method = "complete") 
  
  plot_list[[a]] = x[[4]]     ##to save each plot into a list. note the [[4]]
  
}
g<- grid.arrange(arrangeGrob(grobs= plot_list,ncol=3))
ggsave(paste("output/","RNA-seq_categorized_heatmap_",show_trend,"_v4",".pdf",sep =""),g, width =15,height = 22)