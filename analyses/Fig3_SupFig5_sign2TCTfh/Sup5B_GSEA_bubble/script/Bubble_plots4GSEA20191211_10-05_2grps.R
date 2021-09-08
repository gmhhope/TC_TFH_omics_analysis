setwd("/Volumes/LM_MG_drive_one/Metabolomics_Transcriptome_integration_06042020/Pathway_wise_analysis/TC_Tfh_metabolomics/GSEA_output_july12/Figure3_bubble_plot/")

file_name = "Combn.csv"
df_filt <- read.csv(file_name, header = TRUE)
colnames(df_filt)
df_filt$log10FDR<- -log(df_filt$FDR.q.val,10)
df_filt[df_filt$log10FDR==Inf,"log10FDR"] = 5
df_filt$group <- factor(df_filt$group, levels = c("Tn","Tfh"))
df_filt$NAME = factor(df_filt$NAME, levels=unique(df_filt$NAME[order(df_filt$NES)]), ordered=TRUE)
#fn = factor(f, levels=unique(f[order(a,b,f)]), ordered=TRUE)
colnames(df_filt)
range(df_filt$log10FDR)
#[1] 1.301833 5.000000
range(df_filt$NES)
#-2.316198  2.139029

dim(df_filt)
#[1] 12  6
library(stringr)
library(ggplot2)
 
    pdf(paste("ggplot_bubble_mod",str_extract(file_name, '^[^.]+'),".pdf", sep = "") , width = 5.6, height = 3.6)
    ggplot(df_filt, aes(x=group, y=NAME)) + #reorder(NAME,df_filt$order)
    geom_point(aes(size=log10FDR, colour = NES)) +
    scale_size_continuous(range = c(2, 6)) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", breaks = seq(-3.5,3.5),1) +
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
    