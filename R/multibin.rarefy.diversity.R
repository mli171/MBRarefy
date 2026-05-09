#' Perform Multi-Bin Rarefy for TCR Diversity Across Replicates
#'
#' Executes multi-bin rarefying procedures repeatedly across multiple replicates, for either alpha diversity
#' estimation, using per-sample input files from a directory. This version reads data directly from disk for each
#' replicate to reduce memory footprint. Both sequential and parallel processing are supported.
#'
#' @param InputDataDir Character. Path to the directory containing one file per sample. Each file should
#'   contain tabular data with at least a sequence column and a count column.
#' @param SeqVar Character. Name of the column containing sequence or feature identifiers (used only for beta diversity).
#' @param CountVar Character. Name of the column containing count data (required for both alpha and beta diversity).
#' @param div.measure Character. Either \code{"alpha"} or \code{"beta"} indicating the type of diversity to compute.
#' @param depths Numeric vector. Rarefaction depths to evaluate, e.g., \code{c(1000, 5000, 10000)}.
#' @param methods Character vector. Names of diversity metrics to calculate. For alpha diversity, these might include
#'   \code{"shannon_alpha"}, \code{"gini_simpson_alpha"}, etc.
#' @param nRep Integer. Number of rarefaction replicates to perform. Default is 1.
#' @param parallel Logical. Whether to run replicates in parallel. Default is \code{FALSE}.
#'   Note: parallelization may increase memory usage and is not recommended if memory is limited.
#' @param nCore Integer. Number of CPU cores to use for parallel computation. Required if \code{parallel = TRUE}.
#'
#' @return A named list of length \code{nRep}, each corresponding to a rarefaction replicate.
#'   Each replicate contains the output structure defined by either \code{\link{rarefy.alpha}}.
#'
#' @section Output Structure:
#' The structure depends on \code{div.measure}.
#' \itemize{
#'   \item If \code{"alpha"}: each replicate contains a nested list of per-sample alpha diversity values.
#' }
#'
#' Example for \code{div.measure = "alpha"}:
#' \preformatted{
#' List of nRep
#' ├── Rep1
#' │   ├── Sample1
#' │   │   ├── RarefyTo1000
#' │   │   │   ├── shannon_alpha
#' │   │   │   ├── gini_simpson_alpha
#' |   |   |   ├── ...
#' │   │   ├── RarefyTo5000
#' |   |   |   ├── ...
#' │   │   └── ...
#' │   └── Sample2
#' │   └── ...
#' ├── Rep2
#' └── ...
#' }
#'
#' @details This function wraps around either \code{\link{rarefy.alpha}},
#' based on the \code{div.measure} argument. When \code{parallel = TRUE},
#' it uses \pkg{foreach} and \pkg{doParallel}
#' to distribute replicates across cores. If any required arguments are missing
#' or misconfigured, the function halts with an error.
#'
#' It uses file-based sample-wise processing to reduce cohort-level memory
#' requirements. For alpha diversity, the current rarefaction backend performs
#' exact subsampling without replacement by expanding counts within each sample,
#' so memory use for a single rarefaction task depends on that sample's
#' library size.
#'
#' @seealso \code{\link{rarefy.alpha}}.
#' @importFrom foreach %dopar%
#' @examples
#' \dontrun{
#' set.seed(123)
#' readsLists <- list(
#'   SampleA = list(Seq = c("A", "B", "C"), Rd = c(20, 10, 5)),
#'   SampleB = list(Seq = c("A", "D", "E"), Rd = c(15, 15, 10)),
#'   SampleC = list(Seq = c("C", "D", "F"), Rd = c(5, 10, 15))
#' )
#' depths <- c(1000, 5000, 10000)
#'
#' alphamethods  <- c(
#'   "unique_seq_alpha",      # Number of unique taxa (Richness)
#'   "shannon_alpha",         # Shannon entropy
#'   "gini_simpson_alpha",    # Gini-Simpson index
#'   "chao1_alpha",           # Chao1 richness estimator
#'   "chao1_bc_alpha",        # Bias-corrected Chao1
#'   "pielou_alpha"           # Pielou's evenness
#' )
#'
#' result <- multibin.rarefy.diversity(
#'   InputDataDir = InputDataDir,
#'   SeqVar = "seq",
#'   CountVar = "count",
#'   div.measure = "alpha",
#'   depths = depths,
#'   methods = alphamethods,
#'   nRep = 4,
#'   parallel = TRUE,
#'   nCore = 4
#' )
#' }
#' @export
multibin.rarefy.diversity = function(InputDataDir,
                                     SeqVar=NULL,
                                     CountVar,
                                     div.measure,
                                     depths,
                                     methods,
                                     nRep=1,
                                     parallel=FALSE,
                                     nCore=NULL){

  if(div.measure == "alpha"){
    div.fun = rarefy.alpha
  }else if (div.measure == "beta"){
    # div.fun = rarefy.beta
    stop("Beta diversity is not implemented in this version of MBRarefy.")
  }

  if(parallel){
    nAvaCore = parallel::detectCores()
    if(is.null(nCore)){
      stop(paste0("Missing number of computing cores (",
                  nAvaCore, "cores available)."))
    }
    doParallel::registerDoParallel(cores = nCore)
    divLists = foreach::foreach(i=1:nRep) %dopar%
      (
        do.call(div.fun,
                c(list(InputDataDir, SeqVar, CountVar, depths, methods)))
      )
    names(divLists) = paste0("Rep", 1:nRep)
  }else{
    divLists = vector("list", nRep)
    names(divLists) = paste0("Rep", 1:nRep)
    for(iii in 1:nRep){
      divLists[[iii]] = do.call(div.fun,
                                c(list(InputDataDir, SeqVar, CountVar, depths, methods)))
    }
  }

  return(divLists)
}
