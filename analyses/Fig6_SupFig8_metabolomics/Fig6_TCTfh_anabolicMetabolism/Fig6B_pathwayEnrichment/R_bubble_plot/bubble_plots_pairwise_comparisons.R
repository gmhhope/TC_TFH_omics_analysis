setwd("/Volumes/LM_MG_drive_one/Morel_Metabolomics_comb_1st2nd_spTnTfh_052620/Rm_C2_TC-5/X5_pathway_analyses/R_bubble_plot/")

file_name = "selected_path_enriched.csv"
df_filt <- read.csv(file_name, header = TRUE)
#df_filt <- df_filt[df_filt$FDR < 0.25,]
colnames(df_filt)
df_filt$log10FDR<- -log(df_filt$FDR,10)
#df_filt[df_filt$log10FDR==Inf,"log10FDR"] = 5
levels(df_filt$group)
#"B6_TfhvsB6_Tn"  "TC_TfhvsB6_Tfh" "TC_TfhvsTC_Tn"

df_filt$group <- factor(df_filt$group, levels = c("TC_TfhvsB6_Tfh","TC_TfhvsTC_Tn","B6_TfhvsB6_Tn"))
df_filt = df_filt[order(df_filt$log10FDR,decreasing = TRUE),]
df_filt$Unnamed..0 = factor(df_filt$Unnamed..0, levels=unique(df_filt$Unnamed..0[order(df_filt$log10FDR,decreasing = FALSE)]), ordered=TRUE)
#fn = factor(f, levels=unique(f[order(a,b,f)]), ordered=TRUE)
colnames(df_filt)
range(df_filt$log10FDR)
#[1] 0.6716204 1.8961963

#Calculate fold enrichment
df_filt$Fold_Enrichment = df_filt$hits/df_filt$expected

dim(df_filt)
#[1] 12  6
library(stringr)
library(ggplot2)

pdf(paste("ggplot_bubble_all_red",str_extract(file_name, '^[^.]+'),".pdf", sep = "") , width = 5, height = 3)
ggplot(df_filt, aes(x=group, y=Unnamed..0)) + #reorder(NAME,df_filt$order)
  geom_point(aes(size=Fold_Enrichment, colour = log10FDR)) +
  #scale_size_continuous(range = c(2, 6)) +
  scale_colour_gradient2(low = "#FFCCCB", mid = "red", high = "#8B0000", na.value = NA, midpoint = 1.301,) + #
  #scale_colour_gradient2(low = "blue", mid = "green", high = "red", breaks = seq(0.6,1.3,2)) +
  # set transparency
  # https://ggplot2.tidyverse.org/reference/theme.html
  theme(
    panel.grid.major = element_line(colour = "grey50",linetype = "dashed", size = 0.2),
    panel.border = element_rect(fill = NA),
    #panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    axis.text.x = element_text(angle = 45, hjust = 1,colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
  ) 
dev.off()
