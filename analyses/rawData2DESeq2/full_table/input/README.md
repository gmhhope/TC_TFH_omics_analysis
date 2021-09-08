# RNA-seq full table

# Input
├── README.md
├── result_B6TfhvsB6Tn4grps.csv
├── result_LRT_4grps.csv
├── result_TCTfhvsB6Tfh4grps.csv
├── result_TCTfhvsTCTn4grps.csv
├── result_TCTnvsB6Tn4grps.csv
├── result_interaction4grps.csv
└── rld_4grps.csv

## Notice
- the name of `result_TCTfhvsTCTn4grps.csv` was changed from `Correctedresult_TCTfhvsTCTn4grps.csv` (a minor issue spotted when original result table was generated) to avoid confusion.

## The result tables generated from DESeq2 analyses
- `result_B6TfhvsB6Tn4grps.csv` comparing B6 Tfh versus B6 Tn
- `result_TCTfhvsB6Tfh4grps.csv` comparing TC Tfh versus B6 Tfh
- `result_TCTfhvsTCTn4grps.csv` comparing TC Tfh versus TC Tn
- `result_TCTnvsB6Tn4grps.csv` comparing TC Tn versus B6 Tn
- `result_LRT_4grps.csv` is Likelihood-ratio test on four groups.
- `result_interaction4grps.csv` investigate interative effects between cell types and genotype.