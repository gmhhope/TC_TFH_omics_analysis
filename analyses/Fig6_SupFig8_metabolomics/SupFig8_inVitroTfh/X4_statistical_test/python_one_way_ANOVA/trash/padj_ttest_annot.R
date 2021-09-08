#This script is for ttest4mcg
# setwd("/Volumes/LM_MG_drive_one/Morel_Metabolomics_spTnTfh_01192020/Statistical_test/python_pairwise_t_test/")
# file_name <- "ttest4mcg__PBSvsX400ng.txt"
# df <- read.table(file_name, header = TRUE)
# df[,3] <- p.adjust(df[,3], method = 'BH')
# colnames(df)[3] <- paste("padj_",colnames(df)[3],sep = "")
# 
# write.table(df,paste("padj_",file_name,sep = ""),sep = '\t', row.names = FALSE)

#This part is doing other metabolomics data without intention for mcg.
setwd("/Volumes/LM_MG_drive_one/Morel_Metabolomics_2nd_set_spTnTfh_04202020/rv_TC-5/X3_statistical_test/python_one_way_ANOVA/")
file_name = "ftest4_B6_TFHvsB6_TNvsTC_TFHvsTC_TN.txt"
df <- read.table(file_name, header = TRUE)
#Limit to annotated metabolites to increase singificance
if(T) {
  df <- df[1:99,]  #99
}
df[,4] <- p.adjust(df[,2], method = 'BH', n=99) #99
colnames(df)[4] <- paste("padj_lim_",colnames(df)[3],sep = "")
write.table(df,paste("padj_on_annot_",file_name,sep = ""),sep = '\t', row.names = FALSE)
