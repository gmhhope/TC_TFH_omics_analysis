ctab_loc_dir =  "/Volumes/LM_MG_drive_one/Morel_RNA-seq_UT_071619received/MG_analysis_071619/B6_TC_2DG_project_11262019/DESeq2/Tfh_B6vsTC_r2/"
ctab_name = "count_mtx_Tfh_B6vsTC_r2.csv"
colData_name = "pheno_metadata_Tfh_B6vsTC_r2.txt"
key_name = "Tfh_B6vsTC_r2"
setwd(ctab_loc_dir)
#ctab = read.table(ctab_name,sep = '\t',header = TRUE,row.names = 1)
ctab = read.csv(ctab_name,header = TRUE,row.names = 1)
colData = read.table(colData_name, sep = '\t', header = TRUE)


# Group	SampleID
# SampleID Group
# 1   LspA_1  LspA
# 2   LspA_2  LspA
# 3   LspA_3  LspA
# 4   LspA_4  LspA
# 5    UF1_1   UF1
# 6    UF1_2   UF1
# 7    UF1_3   UF1
# 8    UF1_4   UF1

library(pheatmap)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = ctab , colData = colData, design =~Group)
dds <- dds[rowSums(counts(dds))>1,]
dds$Group <- relevel(dds$Group, ref="B6")
dds <- DESeq(dds, test = "Wald")
res <- results(dds)
write.csv(res,paste("result_",key_name,".csv", sep = ""))

rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

write.csv(assay(rld),file=paste("rld_",key_name,".csv", sep = ""))
write.csv(assay(vsd),file=paste("vsd_",key_name,".csv", sep = ""))



# resA$padj_Rfunc_fdr <- p.adjust(resA$pvalue, method = "fdr")
# write.csv(resA,"UF1vsLspA_actA_DC_0816_Rpadj.csv")
#Since there are a lot of padj = NA, Do padjust using padj function in R

#bM_num = 5
#res <- results(dds_ward, contrast=c("Group", "Crohn","normal"), independentFiltering=FALSE) 
#res$pvalue[res$baseMean < bM_num] <- NA
#res$padj_new <- p.adjust (res$pvalue,method = "BH")

#write.csv(res,paste("baseMean_" ,bM_num, "_CrohnvsNormal.csv", sep = ""))





#sink("mp.df.Th17.CL.Exp.summary.txt")
#summary(res)
#sink()

#### IHW for DESeq ###
#library(IHW)
#resIHW <- results(dds, filterFun=ihw)
#write.table(resIHW, file="result_IHW.txt", sep = '\t',col.names = TRUE)


pdf(paste("MAplot_",key_name,".pdf",sep=""))
plotMA(res, main="DESeq2")
dev.off()
# pdf("MAplot_glyLspA_vs_PBS.pdf")
# plotMA(resB, main="DESeq2", ylim=c(-2,2))
# dev.off()
# pdf("MAplot_nonglyLspAvsPBS.pdf")
# plotMA(resC, main="DESeq2", ylim=c(-2,2))
# dev.off()




#### resMLE using unshrunken dispersion for NB model statistical test ### 
#resMLE <- results(dds, addMLE=TRUE)
#head(resMLE,4)
#pdf("unshrunken_MAplot_df.Th17.CL.Exp.pdf")
#plotMA(resMLE, MLE=TRUE, main="unshrnken LFC", ylim=c(-2,2))
#dev.off()
#write.table(res, file="result.txt", sep = '\t',col.names = TRUE)


library("pheatmap")
#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=FALSE)[1:20] #Decreasing = TRUE/FALSE
#nt <- normTransform(dds) #### log2 normalization method (used for comparison with rld & vsd ###
#log2.norm.counts <- assay(nt)[select,]
#pheatmap(log2.norm.counts, cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, annotation_cols=df)
#pheatmap(assay(rld)[select,],cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, annotation_cols=df)
#pheatmap(assay(vsd)[select,],cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, annotation_cols=df)

##### Distance matrix plot using dist() ####
library(RColorBrewer)
sampleDist <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDist)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf("distance_matrix.pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDist, clustering_distance_cols = sampleDist, col=colors)
dev.off()

##### plotPCA #### (can use my own function "plotPCA.mystyle.R")
nrow(assay(rld)) #20829
pdf("plotPCA_rld_ntop_14917_exp.pdf", width = 10, height = 10)
plotPCA(rld,intgroup="Group", ntop = 20829 ) # if returnData=TRUE, return the data; if not, return the plot)
dev.off()

##### 

plotPCA.mystyle <-  function(object, intgroup="condition", ntop=10000, returnData=FALSE) #just used for returnData = False
{
  
  # calculate the variance for each gene
  library(genefilter)
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  loading <- pca$rotation[,1:3]
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- colData(object)[[intgroup]]
  # assembly the data for the plot
  ########## Here we just use the pcs object passed by the end user ####
  PCvalue <- data.frame(pca$x, group=group, intgroup.df, name=colnames(object))
  
  
  
  if (returnData) {
    result <- list(PCvalue,loading, percentVar)
    return(result)
  }
}

test <- plotPCA.mystyle(rld,intgroup="Group",returnData = TRUE) 

pdf("bar_plot_PCs.pdf")
barplot(test[[3]])
dev.off()
write.csv(test[[1]], "PCvalue.csv")
write.csv(test[[2]], "loadings.csv")

### Find out the outlier and look at the Cook distance summmary(box plot) ###
par(mar=c(8,5,2,2))
pdf("boxplot_cookdist.pdf")
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

### Dispersion plot and fittign alternatives. A usefl diagnostic
pdf("plotDisepEsts.pdf")
plotDispEsts(dds)
dev.off()




