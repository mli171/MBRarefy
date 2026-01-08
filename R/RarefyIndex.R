#' Perform One-Time Rarefaction on Sequence Counts
#'
#' Performs a single rarefaction (random subsampling without replacement) on a set of sequences and their read counts, returning the resulting sequences and their counts after rarefying to a specified depth.
#'
#' @param xSeq A character vector of sequence names (e.g., CDR3s or ASV identifiers), typically corresponding to each taxon or sequence in the original dataset.
#' @param xRds An integer vector of read counts for the corresponding sequences in \code{xSeq}.
#' @param depth An integer specifying the total number of reads to subsample (rarefy) from the original sample. Must be less than or equal to \code{sum(xRds)}.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{RarefySeq}}{A character vector of the sequence names that remained after rarefying.}
#'   \item{\code{RarefyRds}}{An integer vector of read counts corresponding to the sequences in \code{RarefySeq}.}
#' }
#'
#' @details
#' This function simulates one round of rarefaction by treating the input counts as repeated sequence indices,
#' then randomly sampling a specified number of reads. It returns the resulting sequence names and the counts after rarefaction.
#'
#' @examples
#' set.seed(123)
#' xSeq = c("SEQ_A", "SEQ_B", "SEQ_C", "SEQ_D")
#' xRds = c(10, 5, 20, 15)  # total = 50
#' depth = 30
#' result = RarefyIndex(xSeq, xRds, depth)
#' print(result)
#'
#' @export
RarefyIndex = function(xSeq, xRds, depth){

  nX = length(xRds)
  s = sum(xRds)
  tempVec = rep(1:nX, xRds)

  # rarefaction
  tempSample = sample(tempVec, depth)
  tempTable = table(tempSample)

  RarefySeq = xSeq[as.numeric(names(tempTable))]
  RarefyRds = as.vector(tempTable)

  return(list(RarefySeq=RarefySeq, RarefyRds=RarefyRds))
}
