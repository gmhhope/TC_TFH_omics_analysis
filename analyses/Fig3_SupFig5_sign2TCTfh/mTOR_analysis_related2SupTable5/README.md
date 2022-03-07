## A mTOR-focusing analysis
To ask whether there is an overlap between mTOR activator and suppressor genes (Long et al. 2021) and genes that are DE in Tfh cells (adjusted p value < 0.05) between normal B6 mice and lupus-prone TC mice, we used a custom R script to detect significant overlap between gene lists. We found that the 47 gene overlap between mTOR activator genes as identified by Long et al. (2021) and genes upregulated in TC vs. B6 TFH cells was significant (p < 0. 00099). The 37 gene overlap between mTOR suppressor genes (Long et al. 2021) and genes downregulated in TC vs. B6 TFH cells was also significant (p < 0.0099). 

## Methods
We used the genome wide annotation for Mouse available in the R package org.Mm.eg.db to define the number of total number genes in the universe for determining whether the overlap between the gene sets described above was significantly greater than expected under the hypergeometric distribution.




