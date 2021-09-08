# These scripts document steps to have annotation in the matrix used to plot in a heatmap

## 1st step
- `create_mean4each_group.ipynb`, add mean information as columns. In one version, we have tried to plot the data using mean of regularized-log transformed values. But then we think it is still better to use original values to represent the variations within each gene.
  - Input: `../../../../rawData2DESeq2/full_table/full_RNAseq_LRT_added.csv`
  - output: `full_RNAseq_LRT_added_mean_added.csv`

## 2nd step
- `annotation_files4heatmap.ipynb` is a similar script like figure 1A (`../../../../Fig1_TranscriptomeCrossPlatform/preprocessing4Fig1A/scripts/annotation_files4heatmap.ipynb`). The details can read the README.md in that analysis.
- It basically add annotation of statistics, pathway of interest, genes of interest in the expression matrix.