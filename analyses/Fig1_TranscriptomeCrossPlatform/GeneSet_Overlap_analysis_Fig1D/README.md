# GSEAoverlap test
- This is the project file running R dotplot for Gene overlap of FC 1.5 padj < 0.05 in both TC vs B6 (mArray & RNA-seq)

## Folder structure
```
├── GSEAoverlap_dot_plot.Rproj
├── README.md
├── input
│   ├── GL_padj0.05_FC1.5_mArray_RNA-seq_spTfh.txt
│   ├── canonical_pathway_padj0.05_FC1.5.csv
│   └── canonical_pathway_padj0.05_FC1.5.xls
├── output
│   └── dotplot_GSEAoverlap_v3.pdf
└── scripts
    └── GSEAoverlap_dot_plot_R.R
```


## input 
- The input here is using gene list `GL_padj0.05_FC1.5_mArray_RNA-seq_spTfh.txt`, which is computed during Venn Diagram step. It was then search against canonical pathway database in MSigDB (http://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp).
- The output is downloaded as `canonical_pathway_padj0.05_FC1.5.xls` and under further cleanup to produce `canonical_pathway_padj0.05_FC1.5.csv`. 

# workflow
- `GSEAoverlap_dot_plot_R.R` was used to generate the plot and redundant pathways were manually removed to get a succint presentation.