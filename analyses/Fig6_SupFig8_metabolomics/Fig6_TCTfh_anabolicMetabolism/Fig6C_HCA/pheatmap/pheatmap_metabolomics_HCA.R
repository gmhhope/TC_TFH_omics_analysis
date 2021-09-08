
#key option need to modify
#change the number of the pathway, since there are several pathways; run just like manually run loop.
dir = "/Volumes/LM_MG_drive_one/Morel_Metabolomics_comb_1st2nd_spTnTfh_052620/Rm_C2_TC-5/X3_heatmap_HCA_annotated/pheatmap"
setwd(dir)

file_name <- "../m_annot_data_stat_table.csv"
HT <- read.csv(file_name, header = TRUE,row.names = 4)
HT <- HT[HT$padj_lim_ftest_fscore < 0.05,]
colnames(HT) # 11:42 datapoint

colnames(HT)[11:42]
# [1] "C2_B6.1_Tfh_007" "C2_B6.2_Tfh_008" "C2_B6.4_Tfh_009" "C2_B6.6_Tfh_010" "C2_B6.7_Tfh_011" "C2_B6.8_Tfh_012"
# [7] "C2_B6.1_Tn_001"  "C2_B6.2_Tn_002"  "C2_B6.4_Tn_003"  "C2_B6.6_Tn_004"  "C2_B6.7_Tn_005"  "C2_B6.8_Tn_006" 
# [13] "C1_B6.1_Tfh_004" "C1_B6.1_Tn_001"  "C1_B6.2_Tfh_005" "C1_B6.2_Tn_002"  "C1_B6.3_Tfh_006" "C1_B6.3_Tn_003" 
# [19] "C1_TC.1_Tfh_010" "C1_TC.1_Tn_007"  "C1_TC.2_Tfh_011" "C1_TC.2_Tn_008"  "C1_TC.2_Tn_009"  "C1_TC.3_Tfh_012"
# [25] "C2_TC.1_Tfh_018" "C2_TC.2_Tfh_019" "C2_TC.3_Tfh_020" "C2_TC.4_Tfh_021" "C2_TC.1_Tn_013"  "C2_TC.2_Tn_014" 
# [31] "C2_TC.3_Tn_015"  "C2_TC.4_Tn_016" 

colnames(HT)[65] # [1] "padj_lim_ftest_fscore"
retrieve_name <- function(colnames) {
  temp <- sapply(colnames,function(x)(strsplit(x,"_")))
  res <- sapply(temp, function(x)paste(x[3],substr(x[2],1,2), sep = "_"))
  return(res)
}

col_annot <- retrieve_name(colnames(HT)[11:42])
treatment_grp <- unique(col_annot) # [1] "Tfh_B6" "Tn_B6"  "Tfh_TC" "Tn_TC" 

padj_thre = 0.05

library(rlist)
library(pheatmap)
library(RColorBrewer)

zscore <- function(x) {
  z <- (x - rowMeans(x))/ apply(x,1,sd)
  return(z)
}


numeric_input <- HT[,11:42]  #!
z_transformed <-zscore(numeric_input)
z_transformed <- z_transformed[complete.cases(z_transformed), ]

# Add annotation as described above, and change the name of annotation
annotation_col = data.frame(Genotype = col_annot) 
rownames(annotation_col) = colnames(numeric_input)

# annotation_row = data.frame(GeneClass = factor(HT$annotation)) 
# rownames(annotation_row) = rownames(numeric_input)
#Specify column colors
Genotype = c("#0000ff","#000088", "#ff0000","#880000")
names(Genotype) = treatment_grp # [1] "Tfh_B6" "Tn_B6"  "Tfh_TC" "Tn_TC" 
#Output of Genotype
# > Genotype
# Tfh_B6     Tn_B6    Tfh_TC     Tn_TC 
# "#0000ff" "#000088" "#ff0000" "#880000" 

#Specify row colors
# path = list()
# path_colors = colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(length(unique(HT$pathways)))
# path_colors_factor = unique(HT$pathways)
# for (i in (1:length(HT$pathways))) {
#   for (j in (1:length(path_colors_factor))) {
#     if (HT$pathways[[i]]==path_colors_factor[[j]]) {path[[i]] = path_colors[[j]]}
#   }
# }
#names(path) = HT$pathways


#put colors of row and column together
#ann_colors = list(Genotype = Genotype, Path = path)

range(z_transformed) #check the range of z_transformed
breaksList = seq(-2.4, 2.4, 0.1)
cell_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(breaksList)) #rev()


mat_breaks <- seq(min(z_transformed), max(z_transformed), length.out = 1)

out_pdf_file = paste("metabolomics_",".pdf",sep = "")  #gsub(".csv", "",  file_name )

pdf(file=out_pdf_file, width =15,height=20, paper = "special",onefile=FALSE)
# Either z_transformed or HT #
pheatmap(z_transformed, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = annotation_col, clustering_method = "complete", breaks = breaksList, color = cell_colors, cellwidth = 10, cellheight = 10) 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()

















