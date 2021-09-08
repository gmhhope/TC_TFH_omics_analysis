ctab_loc_dir =  "/Volumes/LM_MG_drive_one/Morel_RNA-seq_UT_071619received/MG_analysis_071619/B6_TC_2DG_project_11262019/DESeq2/4grps_TC_B6_Tfh_Tn_interaction//"
ctab_name = "count_mtx_B6_Tfh_vs_Tn_TC_Tfh_vs_Tn.csv"
colData_name = "pheno_metadata_4grp.txt"
key_name = "4grps"
setwd(ctab_loc_dir)
#ctab = read.table(ctab_name,sep = '\t',header = TRUE,row.names = 1)
ctab = read.csv(ctab_name,header = TRUE,row.names = 1)
colData = read.table(colData_name, sep = '\t', header = TRUE)

library(pheatmap)
library(DESeq2)
library(dplyr)
library(tidyr)

#SampleID Genotype CellType
dds<-DESeqDataSetFromMatrix(countData=ctab,colData=colData,design = ~Genotype+CellType+Genotype:CellType);
dds$Genotype <- relevel(dds$Genotype, ref="B6") #It is critical to relevel for interaction terms
dds$CellType <- relevel(dds$CellType, ref="Tn") #It is critical to relevel for interaction terms
dds<-DESeq(dds)
matrix(resultsNames(dds)) #Check options
#     [,1]                    
#[1,] "Intercept"             
#[2,] "Genotype_TC_vs_B6"     
#[3,] "CellType_Tfh_vs_Tn"    
#[4,] "GenotypeTC.CellTypeTfh"



#TC,B6,Tfh,Tn produce result tables with different contrasts
#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
dds$group <- factor(paste0(dds$Genotype, dds$CellType))
design(dds) <- ~ group
dds <- DESeq(dds)
matrix(resultsNames(dds))



dds$group <- relevel(dds$group, ref="B6Tn") #It is critical to relevel for interaction terms
dds <- DESeq(dds)
res2 = results(dds, contrast=c("group", "B6Tfh", "B6Tn"))
write.csv(res2,paste("result_B6TfhvsB6Tn",key_name,".csv", sep = ""))

dds$group <- relevel(dds$group, ref="TCTn") #It is critical to relevel for interaction terms
dds <- DESeq(dds)
res3 = results(dds, contrast=c("group", "TCTfh", "TCTn"))
write.csv(res3,paste("result_TCTfhvsTCTn",key_name,".csv", sep = ""))

dds$group <- relevel(dds$group, ref="B6Tfh") #It is critical to relevel for interaction terms
dds <- DESeq(dds)
res4 = results(dds, contrast=c("group", "TCTfh", "B6Tfh"))
write.csv(res4,paste("result_TCTfhvsB6Tfh",key_name,".csv", sep = ""))

dds$group <- relevel(dds$group, ref="B6Tn") #It is critical to relevel for interaction terms
dds <- DESeq(dds)
res5 = results(dds, contrast=c("group", "TCTn", "B6Tn"))
write.csv(res5,paste("result_TCTnvsB6Tn",key_name,".csv", sep = ""))


#produce result tables with interaction term
dds<-DESeqDataSetFromMatrix(countData=ctab,colData=colData,design = ~Genotype+CellType+Genotype:CellType);
dds$Genotype <- relevel(dds$Genotype, ref="B6") #It is critical to relevel for interaction terms
dds$CellType <- relevel(dds$CellType, ref="Tn") #It is critical to relevel for interaction terms
dds$group <- factor(paste0(dds$Genotype, dds$CellType))
dds<-DESeq(dds)
matrix(resultsNames(dds)) #Check options

res<-results(dds, name="GenotypeTC.CellTypeTfh") #Produce the interaction term
write.csv(res,paste("result_interaction_",key_name,".csv", sep = ""))

res05 <- results(dds, name="GenotypeTC.CellTypeTfh", alpha = 0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
#535

padj_threshold = 0.05
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < padj_threshold)
write.csv(resSig,paste("result_interaction_padj0.05",key_name,".csv", sep = ""))

dds$group <- factor(paste0(dds$Genotype, dds$CellType))

pdf("test_plotCount.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="group")
dev.off()


rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

write.csv(assay(rld),file=paste("rld_",key_name,".csv", sep = ""))
write.csv(assay(vsd),file=paste("vsd_",key_name,".csv", sep = ""))


# plot heatmap Coefficient = log2FC; Unbelievable!
betas <- coef(dds)
colnames(betas)
colnames(betas)
topGenes <- head(order(res$padj),20)
mat <- betas[topGenes, ]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pdf(file="test.pdf", width =5,height=5, paper = "special",onefile=FALSE)
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
dev.off()

#distance matrix
library(RColorBrewer)
sampleDist <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDist)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf("distance_matrix.pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDist, clustering_distance_cols = sampleDist, col=colors)
dev.off()
