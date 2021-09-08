# SupFigure 1F

## directory structure
```
├── Heatmap.Rproj
├── README.md
├── input
│   ├── 3comparisons_AND_cleaned_combn4further_analyses_padj0.1.csv
│   ├── either_R_N_OR_M_N_cleaned_combn4further_analyses_padj0.1.csv
│   └── nanostring_mArray_RNAseq_table062520.xlsx
├── output
│   ├── Heatmap_quantiled_AND_3compare_padj0.1.pdf
│   ├── Heatmap_quantiled_either_Nano_Rseq_or_Nano_mArray_padj0.1.pdf.pdf
│   └── v2_Edited_Heatmap_quantiled_either_Nano_Rseq_or_Nano_mArray_padj0.1.pdf
└── scripts
    └── R-heatmap_on_all_genes_quantiled_coloring_annot_Nanostring_RNAseq_mArray.R
```

## input
- The nanostring data is from Choi et al (Nature Communications 2018).
- The statistical table was merged with RNA-seq and microarray statistical table: `nanostring_mArray_RNAseq_table062520.xlsx`
- The table was then filtered in two ways:
  - DEGs in either (RNA-seq & Nanostring) or (microarray & Nanostring) comparison are significantly (adjusted p-value < 0.1) co-regulated in the same direction between TC Tfh versus B6 Tfh. `either_R_N_OR_M_N_cleaned_combn4further_analyses_padj0.1.csv`
  - We also gathered significant DEGs of which expression in all three dataset share the same trend. `3comparisons_AND_cleaned_combn4further_analyses_padj0.1.csv`

## workflow and output
- The `either_R_N_OR_M_N_cleaned_combn4further_analyses_padj0.1.csv` was used for plotting heatmap in supplementary figure 1F by using the script `R-heatmap_on_all_genes_quantiled_coloring_annot_Nanostring_RNAseq_mArray.R`

