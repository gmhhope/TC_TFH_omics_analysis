library(purrr)
library(tidyverse)
library(dplyr)
library(org.Mm.eg.db)

#The changes from v2 are:
#Only looking at pipeline 2 high stim log2FC from Long et al.
#Looking at overlap in both directions

#Define overlap significance function using hypergeometric distribution
overlap_significance <- function(genes_all, gene_sets, iterations) {
  observed <- length(purrr::reduce(gene_sets, intersect))
  simulated <- purrr::map_dbl(seq_len(iterations), function(x) {
    sim <- purrr::map(lengths(gene_sets), ~sample(genes_all, .x))
    sim <- length(purrr::reduce(sim, intersect))
    return(sim)
  })
  pval <- (sum(simulated >= observed) + 1) / (iterations + 1)
  return(list(pval=pval, simulated_values=simulated, observed=observed))
}

#Read in all mouse genes (symbols)
all_genes <- (keys(org.Mm.eg.db,keytype="SYMBOL"))

##############################
##Reading and filtering data##
##############################
#Read in data
DESeq2=read.csv("Input/full_RNAseq_LRT_added.csv")
#Set cutoffs for overlap significant testing.
log2FCcutoff=0 
padjcutoff=0.05 

#The naming is Inunit.thisunitisdownorup.direction
#UPREG
#TCTfhvsB6Tfh up
Tfh.TC.up=DESeq2 %>% dplyr::select(gene_symbol,log2FoldChange_TCTfhvsB6Tfh,padj_TCTfhvsB6Tfh)%>%
  filter(padj_TCTfhvsB6Tfh<padjcutoff & log2FoldChange_TCTfhvsB6Tfh> log2FCcutoff)
#TCTfhvsB6Tfh down
Tfh.TC.down=DESeq2 %>% dplyr::select(gene_symbol,log2FoldChange_TCTfhvsB6Tfh,padj_TCTfhvsB6Tfh)%>%
  filter(padj_TCTfhvsB6Tfh<padjcutoff & log2FoldChange_TCTfhvsB6Tfh< -log2FCcutoff)

#Read in genes from mTOR paper. #The down and up naming is changed from v1
mTOR.1=read.csv("Input/Long_pipeline1.csv",header=TRUE)
colnames(mTOR.1)=c("gene_symbol","log2FC.pS6highvspS6low.p1.lowstim","log2FC.pS6highvspS6low.p1.highstim")
mTOR.2=read.csv("Input/Long_pipeline2.csv",header=TRUE)
colnames(mTOR.2)=c("gene_symbol","log2FC.pS6highvspS6low.p2.lowstim","log2FC.pS6highvspS6low.p2.highstim")
mTOR=merge(mTOR.1,mTOR.2,by="gene_symbol")

#Pipeline 2
mTOR.p2.highstim.down=mTOR %>% dplyr::select(gene_symbol,log2FC.pS6highvspS6low.p2.highstim) %>% filter(log2FC.pS6highvspS6low.p2.highstim>log2FCcutoff)
mTOR.p2.highstim.up=mTOR %>% dplyr::select(gene_symbol,log2FC.pS6highvspS6low.p2.highstim) %>% filter(log2FC.pS6highvspS6low.p2.highstim< -log2FCcutoff)


#############################
##Running overlap function###
#To remove results rm(list=ls(pattern="*overlap$"))#
#############################

#Are Genes that promote mTORC1 activation and over expressed in TC Tfh as compared B6 Tfh cells?
#Look at the overlap between genes that are upregulatd in TC-Tfh relative to B6-Tfh cells
#and genes that promote mTOR activation as defined by Ps6 hi vs. pS6 lo T-reg

#Overlap between genes that are up-regulated in TC genotype and mTOR activator genes
Tfh.TC.up.mTOR.p2.highstim.overlap=overlap_significance(genes_all = all_genes,gene_sets = c(as.data.frame(Tfh.TC.up$gene_symbol),as.data.frame(mTOR.p2.highstim.up$gene_symbol)),iterations = 1000)
#Overlap between genes that are down-regulated in TC genotype and mTOR supressor genes
Tfh.TC.down.mTOR.p2.highstim.overlap=overlap_significance(genes_all = all_genes,gene_sets = c(as.data.frame(Tfh.TC.down$gene_symbol),as.data.frame(mTOR.p2.highstim.down$gene_symbol)),iterations = 1000)

#Overlap significance#
lists=mget(ls(pattern="*overlap$")) #Right now the individual overlap gene lists are not written to output
#loop over the lists to get p value
overlap.pvalues=as.data.frame(lapply(lists, `[[`, 1))
overlap.pvalues=t(overlap.pvalues)
overlap.pvalues=cbind(rownames(overlap.pvalues), data.frame(overlap.pvalues, row.names=NULL))
colnames(overlap.pvalues)=c("overlap","pval")

#Set cutoffs for generating overlap if different (may want less stringent)
#log2FCcutoff=0 
#padjcutoff=0.05

##THE GENE LISTS BELOW ARE BASED ON NO P-VALUE FILTERING (<1) AND FC FILTERING OF 0.088
##THIS IS TO GET THE DATA FOR THE GENES OF INTEREST TO DR. MOREL
#List of genes in overlap between genes that are up-regulated in TC genotype and mTOR activator genes
overlap.up.list=intersect(Tfh.TC.up$gene_symbol,mTOR.p2.highstim.up$gene_symbol)
#List of genes in overlap between genes that are down-regulated in TC genotype and mTOR supressor genes
overlap.down.list=intersect(Tfh.TC.down$gene_symbol,mTOR.p2.highstim.down$gene_symbol)
#Combine in a way that's like cat|sort|uniq
gene_lists=mget(ls(pattern="*list$"))
genes_in_overlaps <- do.call(c,gene_lists)
genes_in_overlaps <- as.data.frame(unique(genes_in_overlaps)) %>% dplyr::rename("gene_symbol"=`unique(genes_in_overlaps)`)
#Now get the FC from the original data
#Gets the Long et al. data for the overlap genes
up.overlap.genes.FC=merge(genes_in_overlaps,mTOR,all.x=TRUE)
up.overlap.genes.FC=merge(up.overlap.genes.FC,Tfh.TC.up)
write.csv(up.overlap.genes.FC,"mTOR_activator_genes_upreg_in_TCTfhvsB6Tfh.v2.csv",quote=FALSE,row.names = FALSE)
#Down-reg
#Gets the Long et al. data for the overlap genes
down.overlap.genes.FC=merge(genes_in_overlaps,mTOR,all.x=TRUE)
down.overlap.genes.FC=merge(down.overlap.genes.FC,Tfh.TC.down)
write.csv(down.overlap.genes.FC,"mTOR_activator_genes_downreg_in_TCTfhvsB6Tfh.v2.csv",quote=FALSE,row.names = FALSE)
#After writing to output, remove mTOR data columns besides p2 high stim and put both csv files as two tabs in a single excel file.


