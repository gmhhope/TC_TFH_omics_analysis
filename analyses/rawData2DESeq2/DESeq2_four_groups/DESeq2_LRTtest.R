ctab_loc_dir =  "/Volumes/LM_MG_drive_one/Morel_RNA-seq_UT_071619received/MG_analysis_071619/B6_TC_2DG_project_11262019/DESeq2/4grps_TC_B6_Tfh_Tn_interaction//"
ctab_name = "count_mtx_B6_Tfh_vs_Tn_TC_Tfh_vs_Tn.csv"
colData_name = "pheno_metadata_4grp.txt"
key_name = "4grps"
setwd(ctab_loc_dir)
#ctab = read.table(ctab_name,sep = '\t',header = TRUE,row.names = 1)
ctab = read.csv(ctab_name,header = TRUE,row.names = 1)
colData = read.table(colData_name, sep = '\t', header = TRUE)

library(pheatmap)
library(DESeq2)
library(dplyr)
library(tidyr)
# The full model was specified previously with the `design = ~ sampletype`:
# dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
dds_lrt <-DESeqDataSetFromMatrix(countData=ctab,colData=colData,design = ~Genotype+CellType+Genotype:CellType);
dds_lrt$group <- factor(paste0(dds_lrt$Genotype, dds_lrt$CellType))
design(dds_lrt) <- ~ group
matrix(resultsNames(dds_lrt))

# Likelihood ratio test
dds_lrt$group <- relevel(dds_lrt$group, ref="B6Tn") #It is critical to relevel for interaction terms
dds_lrt <- DESeq(dds_lrt, test="LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
write.csv(res_lrt,paste("result_LRT_",key_name,".csv", sep = ""))
