# GSEA bubble plot

## folder structure
├── README.md
├── input
│   ├── README.md
│   ├── gset_report_merged_Sp_Tfh_2grp.csv
│   └── gset_report_merged_microarray_Tfh_2grp.csv
├── output
│   ├── ggplot_bubble_modmeta_df4bubble_mod.pdf
│   ├── meta_df4bubble_diff_filtering_select.csv
│   └── meta_df4bubble_mod.csv
└── scripts
    ├── Bubble_plots4GSEA20191211_10-05_2grps.R
    └── GSEA_clean_up4meta_analysis.ipynb

## workflow
- `GSEA_clean_up4meta_analysis.ipynb` takes `gset_report_merged_Sp_Tfh_2grp.csv` & `gset_report_merged_microarray_Tfh_2grp.csv`. It results in `meta_df4bubble_mod.csv` & `meta_df4bubble_diff_filtering_select.csv`.
  - This notebook is to reduce the size of significant gene list, as plotting all of them is impossible. And redundancy is also prevalent in KEGG gene sets. We decide to plot the most relevant gene sets.
  - As most of the gene sets are downregulated (with FDR p < 0.01), we includes two most prominent upregulated gene sets that are with FDR p < 0.05. 

- `Bubble_plots4GSEA20191211_10-05_2grps.R` takes `meta_df4bubble_diff_filtering_select.csv` and create the plot.