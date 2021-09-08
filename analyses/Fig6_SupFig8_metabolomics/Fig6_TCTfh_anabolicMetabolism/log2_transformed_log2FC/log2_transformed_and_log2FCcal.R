setwd("/Volumes/LM_MG_drive_one/Morel_Metabolomics_comb_1st2nd_spTnTfh_052620/Rm_C2_TC-5/X0_log2_transformed_log2FC/")
file_name = "featab_simp_norm_annot_rm_TC-5.csv"
df <- read.csv(paste("../../raw_data/",file_name, sep = ""), row.names = 1)

#check whether any columns contain NA value.
apply(df, 2, function(x) any(is.na(x)))

log2_df <- log(df,2)
write.csv(log2_df,paste(gsub("\\..*","",file_name),"_log2.csv",sep = ""))

log2_df$TC_Tfh_geoMean <- rowMeans(log2_df[grepl("TC.*_Tfh", colnames(log2_df))])
log2_df$B6_Tfh_geoMean <- rowMeans(log2_df[grepl("B6.*_Tfh", colnames(log2_df))])
log2_df$TC_Tn_geoMean <- rowMeans(log2_df[grepl("TC.*_Tn", colnames(log2_df))])
log2_df$B6_Tn_geoMean <- rowMeans(log2_df[grepl("B6.*_Tn", colnames(log2_df))])

log2_df$log2FC_Tfh_TCvsB6 = log2_df$TC_Tfh_geoMean-log2_df$B6_Tfh_geoMean
log2_df$log2FC_Tn_TCvsB6 = log2_df$TC_Tn_geoMean-log2_df$B6_Tn_geoMean
log2_df$log2FC_TC_TfhvsTn = log2_df$TC_Tfh_geoMean-log2_df$TC_Tn_geoMean
log2_df$log2FC_B6_TfhvsTn = log2_df$B6_Tfh_geoMean-log2_df$B6_Tn_geoMean

write.csv(log2_df,paste(gsub("\\..*","",file_name),"_log2_wt_FC.csv",sep = ""))
