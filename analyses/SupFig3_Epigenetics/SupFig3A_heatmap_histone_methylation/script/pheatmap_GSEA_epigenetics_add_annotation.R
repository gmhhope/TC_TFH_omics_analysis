
#key option need to modify
#change the number of the pathway, since there are several pathways; run just like manually run loop.
dir = "/Users/gongm/Documents/projects/LM_project/08312021_manuscript/Focused_figures4DrLM/Sup_Fig3/replot"
#!
setwd(dir)

treatment_grp = c("Tfh_B6","Tfh_TC")
num_grp = c(4,3)
sample_num = sum(num_grp)
padj_thre = 0.05

library(rlist)
library(pheatmap)
library(RColorBrewer)
# library("org.Mm.eg.db") #BiocManager::install("org.Mm.eg.db")

zscore <- function(x) {
  z <- (x - rowMeans(x))/ apply(x,1,sd)
  return(z)
}


# keytypes(org.Mm.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"      
# [8] "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MGI"         
# [15] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"      
# [22] "SYMBOL"       "UNIGENE"      "UNIPROT"     


file_name <- "GO_HISTONE_METHYLATION_HT4R_mod.csv"
HT <- read.csv(file_name ,header = TRUE,row.names = 2)
HT <- HT[HT$padj < 0.05,]

numeric_input <- HT[,1:sample_num+1]  #!
z_transformed <-zscore(log(numeric_input+1,2)) 
z_transformed <- z_transformed[complete.cases(z_transformed), ]

# Add annotation as described above, and change the name of annotation
annotation_col = data.frame(Genotype = factor(c(rep(treatment_grp[1],num_grp[1]),rep(treatment_grp[2],num_grp[2])))) #set up your treatment (class), and needs to be factor 
rownames(annotation_col) = colnames(numeric_input)  #e.g., PUF1_1, PUF_2, Dpmt1

annotation_row = data.frame(GeneClass = factor(HT$annotation)) 
rownames(annotation_row) = rownames(numeric_input) #e.g., Gene A, Gene B
#Specify column colors
Genotype = c("#0000ff","#008000")
names(Genotype) = c(treatment_grp[1],treatment_grp[2]) #!
#Output of Genotype
###   Genotype
###   PUF1      Cpmt      Dpmt 
###  "#008000" "#0000ff" "#ffff00"

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

out_pdf_file = paste(gsub(".csv", "",  file_name ),".pdf",sep = "")  #change the number

pdf(file=out_pdf_file, width =15,height=20, paper = "special",onefile=FALSE)
# Either z_transformed or HT #
pheatmap(z_transformed, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
         annotation_col = annotation_col, annotation_row = annotation_row, clustering_method = "complete", breaks = breaksList, color = cell_colors, cellwidth = 10, cellheight = 10) 
#HT or z_transformed; col, cluster_rows,cluster_cols; 

#Reset_it
dev.off()



















#significant
for (j in (1:length(new_list_file))) {
file_name <- new_list_file[j]
HT <- read.csv(file_name ,header = TRUE,row.names = 1)
HT = HT[HT$padj < padj_thre,]
if (sum(HT$padj < padj_thre, na.rm = TRUE)) 
{
  HT$GeneName <- mapIds(org.Mm.eg.db,
                        keys=sapply(HT$mouseSymbol,function(x)as.character(x)),
                        column="GENENAME",
                        keytype="SYMBOL",
                        multiVals="first")
}
else {next}


for (i in 1:length(row.names(HT))) {
  row.names(HT)[i] <- paste(row.names(HT)[i], "_", HT$GeneName[i], sep = "")
}


numeric_input <- HT[,1:sample_num+1]  #!
z_transformed <-zscore(log(numeric_input+1,2)) 
z_transformed <- z_transformed[complete.cases(z_transformed), ]

# Add annotation as described above, and change the name of annotation
annotation_col = data.frame(Genotype = factor(c(rep(treatment_grp[1],num_grp[1]),rep(treatment_grp[2],num_grp[2])))) #set up your treatment (class), and needs to be factor 
rownames(annotation_col) = colnames(numeric_input)  #e.g., PUF1_1, PUF_2, Dpmt1

#annotation_row = data.frame(GeneClass = factor(HT$pathways)) # metabolism A...
#rownames(annotation_row) = rownames(numeric_input) #e.g., Gene A, Gene B
#Specify column colors
Genotype = c("#0000ff","#008000")
names(Genotype) = c(treatment_grp[1],treatment_grp[2]) #!
#Output of Genotype
###   Genotype
###   PUF1      Cpmt      Dpmt 
###  "#008000" "#0000ff" "#ffff00"

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

out_pdf_file = paste(gsub(".csv", "",  file_name ),"padj_",as.character(padj_thre),".pdf",sep = "")  #change the number

if (dim(z_transformed)[1] > 2) {
  pdf(file=out_pdf_file, width =15,height=30, paper = "special",onefile=FALSE)
  # Either z_transformed or HT #
  pheatmap(z_transformed, cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE, 
           annotation_col = annotation_col, annotation_row = NULL, clustering_method = "complete", breaks = breaksList, color = cell_colors, cellwidth = 10, cellheight = 10) 
  #HT or z_transformed; col, cluster_rows,cluster_cols; 
  
  #Reset_it
  dev.off()
}
else if (dim(z_transformed)[1] == 1) {
  pdf(file=out_pdf_file, width =15,height=30, paper = "special",onefile=FALSE)
  # Either z_transformed or HT #
  pheatmap(z_transformed, cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = TRUE, 
           annotation_col = annotation_col, annotation_row = NULL, clustering_method = "complete", breaks = breaksList, color = cell_colors, cellwidth = 10, cellheight = 10) 
  #HT or z_transformed; col, cluster_rows,cluster_cols; 
  
  #Reset_it
  dev.off()
}











#### Here is a way to check the duplication problem, which may encounter if unfortinately. ####

#1. check file_name
#> file_name
#[1] "GSE40274_FOXP3_VS_FOXP3_AND_PBX1_TRANSDUCED_ACTIVATED_CD4_TCELL_UP_HT4R.csv"

#2. Then read the file without row.names = 1
#HT <- read.csv(file_name ,header = TRUE)

#3. try to put mouseSymbol in row.names(HT); which will let you know whare those duplicated bugs.
#> View(HT)
#> row.names(HT) = HT$mouseSymbol
#Error in `.rowNamesDF<-`(x, value = value) : 
#  duplicate 'row.names' are not allowed
#In addition: Warning message:
#  non-unique values when setting 'row.names': ‘Ifitm1’, ‘Ifitm2’, ‘Ifitm3’, ‘Ifitm6’, ‘Ifitm7’ 

#4. remove it from the original file