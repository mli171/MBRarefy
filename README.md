# MBRarefy
A multi-bin rarefying method for testing the association between alpha and beta diversities with covariates. This package provides the “multi-bin” rarefying approach that partitions samples into multiple bins according to their library sizes, conducts rarefying within each bin for alpha and beta diversity calculations, and performs meta-analysis across the association study results from each bin.

## Installation
You can install the latest version of MBRarefy from Github:

```{r}
devtools::install_github("mli171/MBRarefy")
```

## Alpha diversity association analysis

A meta-analysis of the association between alpha diversity and a covariate using a multi-bin rarefying approach can be performed via this package. Samples are first binned according to their library sizes. Within each bin, multiple rarefactions are conducted to reduce subsampling variation. Alpha diversity metrics are computed from the rarefied samples and averaged across replicates. Bin-specific association tests are then performed using simple linear regression inference. Finally, results across bins are combined using three meta-analysis strategies: equal weighting (Multi-bin-Equal), sample-size weighting (Multi-bin-SSW), and inverse-variance weighting (Multi-bin-IVW).

### Alpha diversity calculation at different rarefying levels

The output is a nested list with the following structure

```
List of nRep
├── Rep1
│   ├── Sample1
│   │   ├── RarefyTo1000
│   │   │   ├── shannon.alpha
│   │   │   ├── gini.simpson.alpha
|   |   |   ├── ...
│   │   ├── RarefyTo5000
|   |   |   ├── ...
│   │   └── ...
│   └── Sample2
│   └── ...
├── Rep2
└── ...
```

## Reference

[1] Li, Mo, Xing Hua, Shuai Li, Michael C. Wu, and Ni Zhao. "A multi-bin rarefying method for evaluating alpha diversities in TCR sequencing data." Bioinformatics 40, no. 7 (2024): btae431. https://doi.org/10.1093/bioinformatics/btae431.

