#' Compute Alpha Diversity Metrics for a Count Vector
#'
#' Calculates one or more alpha diversity metrics from a given count vector.
#' This function is intended for internal use in rarefaction pipelines such as
#' \code{\link{rarefy.alpha}} and \code{\link{multibin.rarefy.diversity}}.
#'
#' @param xx A numeric vector of counts (e.g., species, OTUs, ASVs) after
#' rarefying.
#' @param methods A character vector specifying which alpha diversity metrics
#' to compute.
#'  Common options include:
#'   \itemize{
#'     \item{\code{"shannon_alpha"} — Shannon entropy}
#'     \item{\code{"gini_simpson_alpha"} — Gini-Simpson index}
#'     \item{\code{"unique_seq_alpha"} — Number of unique taxa}
#'     \item{\code{"chao1_alpha"} — Chao1 richness estimator}
#'     \item{\code{"chao1_bc_alpha"} — Bias-corrected Chao1}
#'     \item{\code{"pielou_alpha"} — Pielou's evenness}
#'   }
#'
#' @return A named list where each element corresponds to a requested alpha
#' diversity metric. The values are numeric or integer scalars depending on the
#' metric.
#'
#' @examples
#' counts = c(3, 5, 2, 10, 15)
#' metrics = c("shannon_alpha", "gini_simpson_alpha", "unique_seq_alpha")
#' alpha.diversity.calc(counts, metrics)
#'
#' @export
alpha.diversity.calc = function(xx, methods){

  nMethods = length(methods)
  alpha.methods = vector("list", nMethods)
  names(alpha.methods) = methods
  for(m in 1:nMethods){
    method = methods[m]
    if(!is.function(method)){
      alpha.calc = get(method)
    }else{
      stop("Require to define Alpha diversity calculation on ", method, "!")
    }
    alpha.methods[[m]] = alpha.calc(xx)
  }

  return(alpha.methods)
}
