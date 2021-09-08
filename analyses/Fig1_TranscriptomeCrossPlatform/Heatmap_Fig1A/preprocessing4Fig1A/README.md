# preoprocessing for Supplementary Figure 1A

## directory structure
├── scripts
│   ├── merge_RNAseq_mArray_result_corrected.html
│   ├── merge_RNAseq_mArray_result.ipynb
│   ├── annotation_files4heatmap.ipynb
│   └── annotation_files4heatmap.html
├── output
│   ├── 2_annotation
│   │   ├── merged_mArray_RNAseq_annot_meta_stat_log2FC_max5_v3.csv
│   │   └── history
│   │       ├── merged_mArray_RNAseq_annot_meta_stat_v2.csv
│   │       └── merged_mArray_RNAseq_annot.csv
│   └── 1_merge_mArray_RNAseq
│       └── merged_mArray_RNAseq.csv
├── input
│   ├── mArray
│   │   ├── Tfh\ tc\ vs\ b6\ t-test\ all\ genes_corrected_01172020.csv
│   │   ├── README.md
│   │   └── Marray_exp_mtx_from_Anton.csv
│   ├── gene_list4annotatition
│   │   ├── gl_KEGG_RIBOSOME_HT4R.csv
│   │   ├── gl_KEGG_OXIDATIVE_PHOSPHORYLATION_4HTprep.csv
│   │   ├── gl_KEGG_CITRATE_CYCLE_TCA_CYCLE_4HTprep.csv
│   │   ├── Tfh_specific_genes_Immunity2014_v2.csv
│   │   └── README.md
│   └── RNA-seq
│       ├── result_Tfh_B6vsTC_r2.csv
│       └── README.md
└── README.md





## input
- README.md explains the detail of the input

## scripts
- `merge_RNAseq_mArray_result.ipynb` is used to merge RNA-seq data with microarray data
- `annotation_files4heatmap.ipynb` is used to annotate the genes and pathway information in the heatmap
- `.html` is their respective exported report

## output
- `1_merge_mArray_RNAseq` is the folder contains the result from `merge_RNAseq_mArray_result.ipynb`
- `2_annotation` is the folder that include the results from `annotation_files4heatmap.ipynb`
  - `merged_mArray_RNAseq_annot_meta_stat_log2FC_max5_v3.csv` 
    - With statistics annotated
      - co-up (no FC cutoff)
      - co-up (at least FC1.5)
      - co-down (no FC cutoff)
      - co-down (at least FC1.5)
    - Restrict to maximum log2FC to 5 in either direction