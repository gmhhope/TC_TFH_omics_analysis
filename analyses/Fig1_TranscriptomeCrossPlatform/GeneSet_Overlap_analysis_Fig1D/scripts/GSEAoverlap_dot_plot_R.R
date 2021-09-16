df = read.csv("./output/combn/combined2dotplot.csv")
df$log10padj = -log(df$FDR.q.value,10)
df$size = as.numeric(df$X..Genes.in.Overlap..k.)
df$Gene.Set.Name = factor(df$Gene.Set.Name, levels = df$Gene.Set.Name)
#df$Gene.Set.Name = factor(df$Gene.Set.Name, levels=unique(sorted_df$Gene.Set.Name[order(df$log10padj,decreasing = FALSE)]),ordered=TRUE)
#data <- read.csv("S1.csv", sep =";", header = TRUE, stringsAsFactors = FALSE)
library("ggplot2")
S1 <- ggplot(df, aes(x=group, y=Gene.Set.Name, size=size, color=log10padj)) + geom_point(alpha = 0.8)  + theme_classic() #+ facet_grid(Group ~ .)
S1

S1 = S1+scale_color_gradient2(low = "white",
                             mid = "gray",
                             high = "red",
                             midpoint = 1.301,
                             space = "Lab",
                             na.value = "grey50",
                             guide = "colourbar",
                             aesthetics = "colour") #limit = c(0.000000000000000000000000000000000000000000000004, 0.03)
S1+scale_size(range = c(2, 8))
pdf("output/GeneSetOverlap_combn_091521.pdf", height = 2.8, width = 6.5)
print(S1)
dev.off()