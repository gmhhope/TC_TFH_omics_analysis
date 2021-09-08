setwd("/Volumes/LM_MG_drive_one/Morel_meta_analyses_transcriptomics/Pearson_correlation_plot_only_TCvsB6_TfhTn_RNAseq/")
df1 = read.csv("result_B6vsTC_r2_Tn.csv")
df2 = read.csv("result_TCvsB6_r2_Tfh.csv")

#filter significant genes
if(T) {
  df1 = df1[df1$padj< 0.05,]
  dim(df1) #~3000
  df2 = df2[df2$padj< 0.05,]
  dim(df2) #~8000
}

library(dplyr)
df1_d = distinct(df1)
df2_d = distinct(df2)
df_m = merge(df1_d,df2_d,by = 'X',how = 'inner')
library("ggpubr")

pdf("sign_genes_stat_pearson.pdf", width = 4,height = 4)
ggscatter(df_m, x = "stat.x", y = "stat.y", alpha = 0.8, size = 0.5,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Wald stat(Tn)", ylab = "Wald stat(Tfh)")
dev.off()

#log2FC
library("ggpubr")
pdf("sign_genes_log2FC_pearson.pdf", width = 4,height = 4)
ggscatter(df_m, x = "log2FoldChange.x", y = "log2FoldChange.y", alpha = 0.8, size = 0.5,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Log2FC_Tn", ylab = "Log2FC_Tfh")
dev.off()