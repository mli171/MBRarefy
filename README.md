# MBRarefy
`MBRarefy` provides an R workflow for alpha diversity association analysis under heterogeneous library sizes in high-throughput count-profile data, including immune-repertoire and microbiome sequencing datasets. The package implements the multi-bin rarefying framework, in which samples are partitioned into library-size bins, rarefied within bins to bin-specific depths, and analyzed by bin-wise association testing followed by cross-bin meta-analysis.

A key feature of `MBRarefy` is automated, data-adaptive library-size cutpoint selection using GA-based ordered knot placement. The package supports fixed-K and varying-K bin selection, repeated rarefying for Monte Carlo stabilization, residual library-size diagnostics, and standardized outputs for downstream association analysis.

## Installation
You can install the latest version of MBRarefy from Github:

```{r}
pak::pak("mli171/MBRarefy")
```

or once available on CRAN

```{r}
install.packages("MBRarefy")
```

## Alpha diversity association analysis

A meta-analysis of the association between alpha diversity and a covariate using a multi-bin rarefying approach can be performed via this package. Samples are first binned according to their library sizes. Within each bin, repeated rarefying are conducted to reduce subsampling variation. Alpha diversity metrics are computed from the rarefied samples and averaged across replicates. Bin-specific association tests are then performed using simple linear regression inference. Finally, results across bins are combined using three meta-analysis strategies: equal weighting (Multi-bin-Equal), sample-size weighting (Multi-bin-SSW), and inverse-variance weighting (Multi-bin-IVW).

## Basic workflow

A typical `MBRarefy` analysis consists of the following steps:

1. Prepare per-sample feature-count files and aligned metadata.
2. Compute repeated rarefaction profiles over a user-defined depth grid using `multibin.rarefy.diversity()`.
3. Aggregate replicate-resolved alpha-diversity results into sample-by-depth matrices using `get_alpha_metric_matrix()`.
4. Select library-size cutpoints using fixed-K or varying-K GA-based optimization.
5. Extract bin-anchored alpha diversity values.
6. Perform bin-wise association testing and cross-bin meta-analysis using `multibin.meta.test.alpha()`.
7. Run a residual library-size diagnostic before interpreting biological or ecological associations.


## Main functions

- `multibin.rarefy.diversity()`: Runs repeated rarefying over a user-specified depth grid and computes alpha diversity metrics for each eligible sample.
- `rarefy.alpha()`: Computes alpha diversity metrics after rarefying for a single sample or count profile.
- `get_alpha_metric_matrix()`: Aggregates replicate-resolved rarefying results into sample-by-depth alpha diversity matrices.
- `fixBinRegObj()`: Objective function for fixed-K GA-based library-size cutpoint selection.
- `varBinRegObj()`: Objective function for varying-K GA-based library-size cutpoint selection.
- `multibin.meta.test.alpha()`: Performs bin-wise alpha diversity association tests and combines bin-specific results by cross-bin meta-analysis.

## Input format

`MBRarefy` uses a file-based input format. Each sample is represented by one plain-text count file with at least a count column and, when required, a feature identifier column. Please check package vignettes how to create input files from commonly seen TCR profiles format and Microbiome OTU count table.

Example:

| seq | count |
|---|---:|
| feature_1 | 15 |
| feature_2 | 3 |
| feature_3 | 27 |

For TCR data, `seq` may represent a CDR3 nucleotide sequence. For microbiome data, `seq` may represent an ASV, OTU, or taxon identifier.

## Quick Start

```{r}
library(MBRarefy)

data(CMV1)

alpha.mats <- get_alpha_metric_matrix(
  CMV1$alpha.res,
  metrics = CMV1$methods
)

USC <- alpha.mats$unique_seq_alpha
head(USC[, 1:3])
```

Calculating alpha diversity at different rarefying levels, the `alpha.mats` object is a nested list with the following structure

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
The complete workflow and demonstrating examples can be found in package vignettes.

The package vignette demonstrates the complete workflow using two application examples:

- a TCR immune-repertoire dataset with known CMV serostatus;
- a wild baboon gut microbiome dataset.

The vignette illustrates repeated rarefying, fixed-\(K\) and varying-\(K\) cutpoint selection, residual library-size diagnostics, and cross-bin meta-analysis.

```{r}
browseVignettes("MBRarefy")
```


## Reference

[1] Li, M. (2026). MBRarefy: data-adaptive multi-bin rarefying for alpha diversity association analysis.

[2] Li, Mo, Xing Hua, Shuai Li, Michael C. Wu, and Ni Zhao. "A multi-bin rarefying method for evaluating alpha diversities in TCR sequencing data." Bioinformatics 40, no. 7 (2024): btae431. https://doi.org/10.1093/bioinformatics/btae431.

