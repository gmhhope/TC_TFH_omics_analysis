setwd("/Volumes/LM_MG_drive_one/2016Immunity_Metabolomics/X2_log2transform/")
file_name = "Jeff-featab_rmTfh5.csv"
df <- read.csv(paste("../X1_Raw_data/",file_name, sep = ""), row.names = 1)

#check whether any columns contain NA value.
apply(df, 2, function(x) any(is.na(x)))

log2_df <- log(df,2)
#scale
scale_log2_df <- t(scale(t(as.matrix(log2_df))))

write.csv(log2_df,paste(gsub("\\..*","",file_name),"_log2.csv",sep = ""))
write.csv(scale_log2_df,paste(gsub("\\..*","",file_name),"_log2_scale.csv",sep = ""))

#log2FC calculation
log2_df$Tact_geoMean <- rowMeans(log2_df[grepl("Tcells", colnames(log2_df))])
log2_df$Tfh_geoMean <- rowMeans(log2_df[grepl("Tfhcell", colnames(log2_df))])
log2_df$Th1_geoMean <- rowMeans(log2_df[grepl("Th1cell", colnames(log2_df))])
log2_df$Th17_geoMean <- rowMeans(log2_df[grepl("Th17cell", colnames(log2_df))])

log2_df$log2FC_TfhvsTact = log2_df$Tfh_geoMean - log2_df$Tact_geoMean
log2_df$log2FC_Th1vsTact = log2_df$Th1_geoMean -log2_df$Tact_geoMean
log2_df$log2FC_Th17vsTact = log2_df$Th17_geoMean - log2_df$Tact_geoMean

write.csv(log2_df,paste(gsub("\\..*","",file_name),"_log2_wt_FC.csv",sep = ""))

