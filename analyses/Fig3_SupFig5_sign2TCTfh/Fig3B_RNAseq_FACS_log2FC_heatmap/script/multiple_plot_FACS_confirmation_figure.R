library(dplyr)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(rlist)
library(gridExtra)
library(ggplot2)
library("org.Mm.eg.db") #BiocManager::install("org.Mm.eg.db")

HT <- read.csv("input/FACS_table.csv", header = T, row.names = 1,check.names = F)
gene_df <- read.csv("input/FACS_RNAseq_gene_name_conv.csv",header = T, stringsAsFactors = TRUE)
category_list <- unique(gene_df$Category) 

colnames(HT)[c(9:12)]# check whether the data columns in this area
data = HT[,c(9:12)]


#annotation of the columns
treatment_grp = c("log2FC_NonTfh_TCvsB6","log2FC_Tfh_TCvsB6","log2FC_B6_TfhvsNonTfh","log2FC_TC_TfhsNonTfh")
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
  scaled = data.matrix(data[gene_df$Category==a,]) #
  
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
g<- grid.arrange(arrangeGrob(grobs= plot_list,ncol=1))
ggsave("output/FACS_categorized_heatmap.pdf",g, width =10,height = 15)