#' Calculate Number of Unique Sequences
#'
#' Computes the number of observed unique sequences, i.e., how many elements
#' in the count vector have nonzero abundance.
#'
#' @param x A numeric vector of counts (e.g., sequence or OTU abundances).
#'
#' @return An integer indicating the number of unique (nonzero) sequences.
#'
#' @examples
#' unique_seq_alpha(c(1, 0, 2, 3, 0))
#' # Returns 3
#'
#' unique_seq_alpha(c(0, 0, 0))
#' # Returns 0
#'
#' @export
unique_seq_alpha = function(x) {
  sum(x > 0)
}



#' Calculate Shannon Diversity Index
#'
#' Computes the Shannon entropy, a commonly used diversity measure in ecology.
#'
#' @param x A numeric vector of counts.
#' @return A numeric value representing the Shannon diversity index.
#' @examples
#' shannon_alpha(c(10, 20, 30))
#' @export
shannon_alpha <- function(x){
  p = x/sum(x)
  return(-sum(p * log(p)))
}



#' Calculate Gini-Simpson Diversity Index
#'
#' Computes the Gini-Simpson index, a measure of diversity representing the probability that two randomly selected individuals belong to different species.
#'
#' @param x A numeric vector of counts.
#' @return A numeric value representing the Gini-Simpson index.
#' @examples
#' gini_simpson_alpha(c(5, 10, 15))
#' @export
gini_simpson_alpha <- function(x){
  p = x/sum(x)
  return(1 - sum(p^2))
}



#' Calculate Chao1 Richness Estimator
#'
#' Computes the Chao1 estimator for species richness based on observed singleton and doubleton species.
#'
#' @param x A numeric vector of counts.
#' @return A numeric value representing the Chao1 estimate of richness.
#' @examples
#' chao1_alpha(c(1, 1, 2, 3, 4))
#' @export
chao1_alpha <- function(x){
  tmpf1 = sum(1*(x == 1))
  tmpf2 = sum(1*(x == 2))
  return(length(x) + (tmpf1^2)/(2*tmpf2))
}



#' Calculate Bias-Corrected Chao1 Richness Estimator
#'
#' Computes a bias-corrected version of the Chao1 richness estimator to account for small sample sizes.
#'
#' @param x A numeric vector of counts.
#' @return A numeric value representing the bias-corrected Chao1 estimate.
#' @examples
#' chao1_bc_alpha(c(1, 1, 2, 3, 4))
#' @export
chao1_bc_alpha <- function(x){
  tmpf1 = sum(1*(x == 1))
  tmpf2 = sum(1*(x == 2))
  return(length(x) + (tmpf1*(tmpf1 - 1))/(2*(tmpf2+1)))
}

#' Calculate Pielou's Evenness Index
#'
#' Computes Pielou's evenness, which quantifies how evenly individuals are distributed across species.
#'
#' @param x A numeric vector of counts.
#' @return A numeric value representing Pielou's evenness index.
#' @examples
#' pielou_alpha(c(10, 10, 10))
#' @export
pielou_alpha <- function(x){
  p = x/sum(x)
  return(-sum(p*log(p))/log(length(x)))
}
