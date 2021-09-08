setwd("/Volumes/LM_MG_drive_one/Morel_Metabolomics_comb_1st2nd_spTnTfh_052620/Rm_C2_TC-5/X4_box_plot/09062021_test_all/")

# library
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

df <- read.csv("../../X2_statistical_test/m_annot_data_stat_table.csv")
metab_list <- read.csv("./metabolite_list_two_way_anova_0.05.csv")
df <- merge(df,metab_list,by.x = "annotated_name_KEGG",by.y = "KEGG_name",how = "inner")
if(FALSE) { # Manual ordering of the plot.
  # order = c("aspartic acid","glutamic acid")
} else {
  order = metab_list$KEGG_name
}


colnames(df)
#df_filt = df[df$padj_lim_ftest_fscore > 0.05 & df$ttest_pval_TC_TfhvsB6_Tfh < 0.05,]
#df_filt = df[df$padj_lim_ftest_fscore > 0.05 & df$ttest_pval_TC_TnvsB6_Tn < 0.05,]
#df <- df_filt

df_sel = df[,c(12:43)] #df_filt
dim(df_sel)
df_t = t(df_sel)
df_t = as.data.frame(df_t)
colnames(df_t) = df$annotated_name_KEGG
length(colnames(df_t))
df_t$group = sapply(strsplit(colnames(df_sel),"_"),function(x)paste(strsplit(x[2],"\\.")[[1]][1],x[3],sep = "_"))  #"." should be either "[.]" or "\\."
#df_t$cohort = sapply(strsplit(colnames(df_sel),"_"),function(x)x[1])

df_long = gather(df_t, metabolites, measurement,1:length(order), factor_key =TRUE) #1:length

df_mean = group_by(df_long, metabolites) %>% summarize(mean = mean(measurement))
df_long = merge(df_long,df_mean, by = "metabolites")
df_long2 = df_long[order(df_long$mean, decreasing = TRUE),]
df_long2$metabolites = factor(df_long2$metabolites,levels  = order) #unique(df_long2$metabolites)
df_long2$group = factor(df_long2$group)
df_long2$group = factor(df_long2$group, levels = c("B6_Tn","TC_Tn","B6_Tfh","TC_Tfh")) #levels = c("B6_Tn","B6_Tfh","TC_Tn","TC_Tfh")
#df_long2$cohort = factor(df_long2$cohort) #factorized

df_long2 = df_long2[order(df_long2$metabolites),]

metab_l = unique(df_long2$metabolites)

#for loop
my_comparisons = list(c("B6_Tfh","TC_Tfh"),c("B6_Tn","TC_Tn"),c("B6_Tfh","B6_Tn"),c("TC_Tfh","TC_Tn"))

fig_l = list()
for (i in 1:length(metab_l)) {
  temp = df_long2[df_long2$metabolites == metab_l[[i]],]
  fig_l[[i]] = ggplot(temp, aes_string(x = "group", y = "measurement")) + 
    geom_boxplot(color=c("black","red","black","red"), fill=(c("#ffffff","#ffffff","black","red")), alpha=0.2) +
    labs(title = metab_l[[i]],
         x ="Groups", y = "log2(intensity)") +
    geom_jitter() + # aes(colour = cohort)Exclude cohort distinguishing criteria
    theme(axis.text.x = element_text(angle = 90, hjust = 1), )+ 
    theme_bw() +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = 't.test', p.adjust.method="fdr")+ #label = "p.signif" ,
    stat_compare_means(method = "anova", label.y = temp$mean[1]-2)     # Add global p-value , 
  
}
multi.page <- ggarrange(plotlist = fig_l,
                        nrow = 3, ncol = 3)
ggexport(multi.page, filename = "two_way_anova_0.05_metabolite_boxplot.pdf", width = 10, height = 8)


# library(ggpubr)
# 
# my_comparisons = list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
# 
# ggboxplot(ToothGrowth, x = "dose", y = "len",
#           color = "dose", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+
#   stat_compare_means(label.y = 45)