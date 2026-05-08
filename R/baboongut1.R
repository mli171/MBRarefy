#' Baboon gut microbiome example dataset
#'
#' A preprocessed example dataset for demonstrating the MBRarefy microbiome
#' alpha-diversity workflow. The dataset is derived from the longitudinal wild
#' baboon gut 16S rRNA microbiome dataset used in the GrieneisenTS study and
#' available through the Bioconductor `microbiomeDataSets` package. Because the
#' original dataset contains repeated longitudinal samples, this example retains
#' one profile per baboon by selecting the sample with the largest library size
#' for each individual, yielding 585 independent baboon samples.
#'
#' The object is used in the package vignette to illustrate grid-based repeated
#' rarefaction, construction of alpha-diversity matrices, genetic-algorithm-based
#' selection of library-size cutpoints, and multibin meta-analysis of associations
#' between observed richness and host or environmental covariates such as age,
#' sex, and season.
#'
#' @format A named list with four elements:
#' \describe{
#'   \item{dataPheno}{A data frame containing sample metadata aligned to the
#'   rarefaction results. Rows correspond to selected baboon samples. Variables
#'   include sample identifiers, baboon identifiers, total read counts
#'   (`totalReads`), and covariates used in the vignette such as age, sex, and
#'   season, when available.}
#'   \item{alpha.res}{A replicate-resolved nested list returned by
#'   \code{multibin.rarefy.diversity()}, containing alpha-diversity estimates
#'   computed over a grid of rarefying depths.}
#'   \item{depths}{Numeric vector of rarefying depths used for grid-based
#'   alpha-diversity profiling.}
#'   \item{methods}{Character vector of alpha-diversity metrics computed, such
#'   as observed richness, Shannon entropy, Gini-Simpson index, Chao1, and
#'   Pielou's evenness.}
#' }
#'
#' @source Derived from the wild baboon gut microbiome dataset available through
#' the Bioconductor `microbiomeDataSets` package. The package vignette describes
#' the subject-level preprocessing and MBRarefy analysis workflow.
#'
#' @references
#' Grieneisen, L. E. et al. (2021). Gut microbiome heritability is nearly
#' universal but environmentally contingent. Science, 373(6551), 181--186.
#'
#' Lahti, L., Ernst, F., and Shetty, S. (2025). microbiomeDataSets:
#' ExperimentHub-based microbiome datasets. Bioconductor package.
#'
#' @examples
#' data(baboongut1)
#' names(baboongut1)
#'
#' dataPheno <- baboongut1$dataPheno
#' depths <- baboongut1$depths
#' methods <- baboongut1$methods
#'
#' head(dataPheno)
#' depths
#' methods
"baboongut1"
