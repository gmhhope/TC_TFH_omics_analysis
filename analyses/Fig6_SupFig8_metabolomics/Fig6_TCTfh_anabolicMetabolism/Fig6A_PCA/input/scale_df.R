setwd("/Volumes/LM_MG_drive_one/Morel_Metabolomics_comb_1st2nd_spTnTfh_052620/Rm_C2_TC-5/X1_PCA/With_scale/")
file_name = "featab_simp_norm_annot_rm_TC-5_log2.csv"
log2_df <- read.csv(paste("./",file_name, sep = ""), row.names = 1)

#scale
scale_log2_df <- t(scale(t(as.matrix(log2_df))))

write.csv(scale_log2_df,paste(gsub("\\..*","",file_name),"_log2_scale.csv",sep = ""))


