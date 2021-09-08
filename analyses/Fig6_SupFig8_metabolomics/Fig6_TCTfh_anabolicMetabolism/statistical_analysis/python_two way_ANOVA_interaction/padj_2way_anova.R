#This script is for ttest4mcg
# setwd("/Volumes/LM_MG_drive_one/Morel_Metabolomics_spTnTfh_01192020/Statistical_test/python_pairwise_t_test/")
# file_name <- "ttest4mcg__PBSvsX400ng.txt"
# df <- read.table(file_name, header = TRUE)
# df[,3] <- p.adjust(df[,3], method = 'BH')
# colnames(df)[3] <- paste("padj_",colnames(df)[3],sep = "")
# 
# write.table(df,paste("padj_",file_name,sep = ""),sep = '\t', row.names = FALSE)

#This part is doing other metabolomics data without intention for mcg.
setwd("/Volumes/LM_MG_drive_one/Morel_Metabolomics_comb_1st2nd_spTnTfh_052620/Rm_C2_TC-5/X2_statistical_test/python_two way_ANOVA_interaction/")
file_name = "2wayAnova.csv"
#df <- read.table(file_name, header = TRUE)
df <- read.csv(file_name, header = TRUE)
#Limit to annotated metabolites to increase singificance
dim(df)
#[1] 99 24
df[,3] <- p.adjust(df[,2], method = 'BH') #
colnames(df)[3] <- paste("padj_lim_",colnames(df)[2],sep = "") #
write.csv(df,paste("padj_on_annot_",file_name,sep = ""),row.names = FALSE)
