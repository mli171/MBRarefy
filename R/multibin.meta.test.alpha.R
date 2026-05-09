#' Multi-bin Meta-analysis for Testing Association Between Alpha Diversity and Covariates
#'
#' Performs a meta-analysis of associations between alpha diversity and a covariate
#' using a multi-bin rarefying approach. Samples are binned by library size, and
#' association tests are performed within each bin. Three meta-analysis weighting
#' strategies are used to combine single-bin results: equal weighting (Equal),
#' sample-size weighting (SSW), and inverse-variance weighting (IVW).
#'
#' @param y.alpha.mat A matrix of alpha diversity values with rows representing
#' samples and columns representing bins.
#' @param xs A numeric vector of covariate values (e.g., binary or continuous) of the same length as the number of samples.
#' @param totalReads A numeric vector of library sizes (e.g., total read counts per sample).
#' @param BinCuts A numeric vector specifying the cut points to bin \code{totalReads}.
#'   Must have length \code{nBins + 1}.
#' @param test.func A function or character name of a function used to test the association
#'   within each bin. The function should accept arguments \code{X} (covariate) and \code{Y} (response).
#'   Defaults to \code{"test.func.bin"}.
#'
#' @return A list with four components:
#' \describe{
#'   \item{\code{Single.bin.res}}{A list of length equal to the number of bins, each containing sample size, test statistic, and p-value from the bin-specific test.}
#'   \item{\code{Multi.bin.Equal}}{Meta-analysis results using equal weights across bins.}
#'   \item{\code{Multi.bin.SSW}}{Meta-analysis using sample-size weights.}
#'   \item{\code{Multi.bin.IVW}}{Meta-analysis using inverse-variance weights.}
#' }
#'
#' @details Each sample is assigned to a bin based on its library size using \code{BinCuts}.
#' Within each bin, a statistical test is applied to evaluate the association between alpha
#' diversity and the covariate. The results are then aggregated using three meta-analytic strategies:
#' \itemize{
#'   \item \strong{Multi-bin-Equal:} Equal weighting of bin-specific estimates.
#'   \item \strong{Multi-bin-SSW:} Sample-size weighted average of bin-specific estimates.
#'   \item \strong{Multi-bin-IVW:} Inverse-variance weighted average of bin-specific estimates.
#' }
#'
#' @references
#' Mo Li, Xing Hua, Shuai Li, Michael C. Wu, and Ni Zhao (2024).
#' A multi-bin rarefying method for evaluating alpha diversities in TCR sequencing data.
#' \emph{Bioinformatics}, 40(7), btae431.
#' \doi{10.1093/bioinformatics/btae431}
#' @importFrom stats na.omit
#' @examples
#' set.seed(123)
#' totalReads = sample(10000:30000, 100, replace = TRUE)
#' xs = rbinom(100, 1, 0.5)
#' BinCuts = seq(10000, 30000, by = 5000)
#' y.alpha.mat = matrix(rnorm(100 * (length(BinCuts) - 1), mean = 7), nrow = 100)
#'
#' # Run multi-bin meta-analysis
#' result = multibin.meta.test.alpha(y.alpha.mat, xs, totalReads, BinCuts)
#'
#' # Access single-bin and IVW results
#' result$Single.bin.res
#' result$Multi.bin.IVW
#'
#' @export
multibin.meta.test.alpha = function(y.alpha.mat,
                                    xs,
                                    totalReads,
                                    BinCuts,
                                    test.func="test.func.bin"){

  nBins = length(BinCuts) - 1

  if(any(totalReads < min(BinCuts) | totalReads > max(BinCuts))){
    stop("Sample library size out of provided Bin range.")
  }

  binFactor = cut(totalReads, BinCuts)

  if (is.function(test.func)) {
    test.bin <- test.func
  } else {
    test.bin <- get(test.func, mode = "function")
  }

  # Single-bin
  Single.bin.res = vector("list", nBins)
  names(Single.bin.res) = paste0("Bin", 1:nBins)
  bClc = vClc = BinSizes = rep(NA, nBins)
  for(i in 1:nBins){
    tempIdxBin = which(!is.na(binFactor) &
                         totalReads >= BinCuts[i] &
                         totalReads <= BinCuts[i+1])
    BinSizes[i] = length(tempIdxBin)
    tmp.bin.test = test.bin(X=xs[tempIdxBin], Y=y.alpha.mat[tempIdxBin,i])
    bClc[i] = as.numeric(tmp.bin.test[1])
    vClc[i] = as.numeric(tmp.bin.test[2]^2)
    Single.bin.res[[i]] = list(Samplesize=length(tempIdxBin),
                               Test.Stat=as.numeric(tmp.bin.test[3]),
                               P.val=as.numeric(tmp.bin.test[4]))
  }

  # Multi-bin-Equal
  wts = rep(1, nBins)
  T.param = mean(bClc)
  T.test.stat = sqrt((sum(bClc))^2 / sum(vClc))
  T.pval = 2 * pnorm(T.test.stat, lower.tail = FALSE)
  Multi.bin.Equal = list(wts=wts,
                         T.param=T.param,
                         T.test.stat=T.test.stat,
                         T.pval=T.pval)

  # Multi-bin-SSW
  wts = BinSizes
  T.param = sum(wts*bClc)/sum(wts)
  T.test.stat = sqrt((sum(wts*bClc)^2) / sum(wts*wts*vClc))
  T.pval = 2 * pnorm(T.test.stat, lower.tail = FALSE)
  Multi.bin.SSW = list(wts=wts,
                       T.param=T.param,
                       T.test.stat=T.test.stat,
                       T.pval=T.pval)

  # Multi-bin-IVW
  wts = (1/vClc)
  T.param = sum(wts*bClc)/sum(wts)
  T.test.stat = sqrt((sum(wts*bClc))^2/ sum(wts*wts*vClc))
  T.pval = 2 * pnorm(T.test.stat, lower.tail = FALSE)
  Multi.bin.IVW = list(wts=wts,
                       T.param=T.param,
                       T.test.stat=T.test.stat,
                       T.pval=T.pval)

  res = list(Single.bin.res = Single.bin.res,
             Multi.bin.Equal=Multi.bin.Equal,
             Multi.bin.SSW=Multi.bin.SSW,
             Multi.bin.IVW=Multi.bin.IVW)

  return(res)
}
