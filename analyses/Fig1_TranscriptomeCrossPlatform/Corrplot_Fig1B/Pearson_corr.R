setwd("/Volumes/LM_MG_drive_one/Morel_meta_analyses_transcriptomics/Pearson_correlation_plot_only_mArray_RNAseq/")
df1 = read.csv("Microarray_res_Tfh_TCvsB6_01172020corrected.csv")
df2 = read.csv("result_TCvsB6_r2_Tfh.csv")

#filter significant genes
if(F) {
  df1 = df1[df1$p.value< 0.01,]
  dim(df1) #~3000
  df2 = df2[df2$padj< 0.05,]
  dim(df2) #~8000
}

library(dplyr)
df1_d = distinct(df1)
df2_d = distinct(df2)
df_m = merge(df1_d,df2_d,by = 'X',how = 'inner')
library("ggpubr")

#check column names
colnames(df_m)

pdf("all_genes_stat_pearson.pdf", width = 4,height = 4)
ggscatter(df_m, x = "t_score_TCvsB6", y = "stat", alpha = 0.8, size = 0.5,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Wald stat(mArray)", ylab = "Wald stat(RNA-seq)")
dev.off()

#log2FC
library("ggpubr")
pdf("all_genes_pearson_log2FC.pdf", width = 4,height = 4)
ggscatter(df_m, x = "log2FC_TCvsB6", y = "log2FoldChange", alpha = 0.8, size = 0.5,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Log2FC_Microarray", ylab = "Log2FC_RNA-seq")
dev.off()