# Supplementary Figure 2F-G

## directory structure
├── R-heatmap_on_selected_genes_quantiled_coloring_Tfhrelated_TfhvsTn_nodiff_btw_TCandB6.R
├── R-heatmap_on_selected_genes_quantiled_coloring_Tfhrelated_TfhvsTn_only_diff_btw_TCandB6.R
├── README.md
├── input
│   └── Tfh_specific_genes_Immunity2014_v3.csv
└── output
    ├── Heatmap_quantiled_01302020_selected_Tfhrelated_TfhvsTn_only_diff_btw_TCandB6.pdf
    └── Heatmap_quantiled_031920_selected_Tfhrelated_TfhvsTn_nodiff_btw_TCandB6.pdf

## input
- The input contains the full table or regularized-logarithm transformed count matrix. In the script, we used the `../../rawData2DESeq2/full_table/full_RNAseq_LRT_added.csv` and identify Tfh related genes in `Tfh_specific_genes_Immunity2014_v3.csv`.
- We sort them out into two part: 
  - one part of genes show no differences between TC and B6 (`Tfhrelated_TfhvsTn_nodiff_btw_TCandB6`)
  - one part of genes of which its expression show significant difference between TC and B^ (`Tfhrelated_TfhvsTn_only_diff_btw_TCandB6`)
  - And their respective Rscript and outputs are in the directory.