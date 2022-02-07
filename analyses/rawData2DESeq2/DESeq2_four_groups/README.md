# DESeq2 output from feature count matrix
- Afte initial pairwise comparisons, I performed the following analyses to facility four-groups comparisons, plotting four groups using the regularized transformed values from DESeq2 and also perform analyses with interaction term.
  - Of note, the initial analyses were in the `../DESeq2_pairwise_comparison/`.
- `DESeq2_LRTtest.R` take cares of likelihood-ratio test among four groups
- `DESeq2_4grp_pairwiseCompare_Interaction_01092020.R` take cares of pairwise comparisons between any two groups and also perform analyses of interaction term between genotype and cell type.
