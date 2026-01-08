#' Compute Alpha Diversity via Rarefaction from Input Files
#'
#' Reads per-sample count tables from an input directory and performs rarefaction-based alpha diversity estimation
#' across multiple depths and metrics. Each sample is independently rarefied at each depth, and alpha diversity is
#' calculated using user-specified metrics.
#'
#' @param InputDataDir Character string specifying the directory path containing input files.
#'   Each file should be a tab-delimited or space-separated count table for a single sample (.txt) file.
#'   The count of each TCR sequences is required for the alpha diversity association analysis.
#' @param SeqVar NULL here. Name of the column containing sequence or feature identifiers (used only for beta diversity
#'   for pairwise comparison).
#' @param CountVar Character string specifying the column name in the input files that contains count data
#'   (e.g., read or taxon counts).
#' @param depths Numeric vector specifying the rarefaction depths to apply (e.g., \code{c(1000, 5000, 10000)}).
#' @param methods Character vector of alpha diversity metric function names to compute for each rarefied sample.
#'   Each method should correspond to a function accepted by \code{\link{alpha.diversity.calc}}.
#'
#' @return A named list of length equal to the number of samples. Each element is itself a list indexed by rarefaction depth
#' (e.g., \code{"RarefyTo1000"}), containing a named list of alpha diversity metrics.
#'
#' @section Output Structure:
#' \preformatted{
#' List of Samples
#' ├── Sample1
#' │   ├── RarefyTo1000
#' │   │   ├── shannon_alpha     : numeric
#' │   │   ├── gini_simpson_alpha: numeric
#' |   |   ├── ...
#' │   ├── RarefyTo5000
#' │   └── ...
#' ├── Sample2
#' └── ...
#' }
#'
#' @details
#' The function reads sample files from the specified directory, extracts the count column, and performs rarefaction
#' at each specified depth (if the sample's total count is at least that depth). It uses multinomial subsampling
#' to select a fixed number of observations and computes diversity metrics from the resulting subsample.
#'
#' If the rarefaction depth exceeds the library size of a sample, the returned metrics for that depth will be set to
#' \code{NULL} (but still structured in the output).
#'
#' This function is commonly used with \code{\link{multibin.rarefy.diversity}} to perform replicate-based evaluation.
#'
#' @seealso \code{\link{alpha.diversity.calc}}, \code{\link{multibin.rarefy.diversity}}
#' @importFrom utils read.csv
#' @examples
#' \dontrun{
#' # Assuming InputDataDir contains one count file per sample
#' InputDataDir <- "data/samples/"
#' CountVar <- "Count"
#' depths <- c(1000, 5000, 10000)
#' methods <- c("shannon_alpha", "gini_simpson_alpha")
#'
#' result <- rarefy.alpha(InputDataDir, CountVar, depths, methods)
#' }
#'
#' @export
rarefy.alpha = function(InputDataDir, SeqVar=NULL, CountVar, depths, methods){

  nBins = length(depths)

  fileNames = list.files(InputDataDir)
  SampleNames = unlist(strsplit(fileNames, ".txt"))
  nSamples = length(SampleNames)

  alpha.samples = vector("list", nSamples)
  names(alpha.samples) = SampleNames
  for(i in 1:nSamples) {

    tempDataRaw = read.csv(file=paste0(InputDataDir, fileNames[i]), sep="")
    ctvec = tempDataRaw[,CountVar]
    rm(tempDataRaw)

    lib.size = sum(ctvec)
    tempvec = rep(1:length(ctvec), ctvec)

    alpha.depths = vector("list", nBins)
    names(alpha.depths) = paste0("RarefyTo", depths)
    for(d in 1:nBins){
      depth = depths[d]
      if(depth > lib.size){
        # NULL value due to library size less than the rarefy level
        alpha.methods = vector("list", length(methods))
        names(alpha.methods) = methods
        alpha.depths[[d]] = alpha.methods
      }else{
        tempSample = sample(tempvec, depth)
        tempTable = table(tempSample)
        if(sum(tempTable) != depth) stop("Rarefying incorrect for Sample: ",
                                         SampleNames[i],
                                         " on Rarefying level: ",
                                         depth, "!")
        alpha.depths[[d]] = alpha.diversity.calc(xx=as.vector(tempTable),
                                                 methods)
      }
    }

    alpha.samples[[i]] = alpha.depths
  }

  return(alpha.samples)
}
