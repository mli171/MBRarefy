#' Compute Summary Matrices for Alpha Diversity Metrics
#'
#' Aggregates alpha diversity metrics across replicates for each sample and
#' rarefaction depth level. The function returns one matrix per metric, with
#' rows corresponding to sample names and columns to rarefaction depths. Each
#' cell represents the average metric value across all replicates.
#'
#' @param x A nested list structure containing alpha diversity values. The
#' expected structure is: `x[[replicate]][[sample]][[depth]][[metric]]`,
#' where `replicate`, `sample`, `depth`, and `metric`
#' are character keys indexing the hierarchy.
#' @param metrics A character vector of metric names to extract and average
#' (e.g., `"shannon_alpha"`, `"unique_seq_alpha"`, `"gini_simpson_alpha"`, etc.).
#'
#' @return A named list of data frames, each corresponding to a metric.
#'   Rows represent samples and columns represent rarefaction depths. Entries
#'   are averages across replicates, with `NA` used for missing or unavailable
#'   values.
#'
#' @examples
#' # Simulated example for alpha diversity summary
#' set.seed(123)
#' # Fake rarefaction results for 2 replicates, 3 samples, 3 depths, 3 metrics
#' reps <- paste0("Rep", 1:2)
#' samples <- c("SampleA", "SampleB", "SampleC")
#' depths <- c("Rarefy1000", "Rarefy5000", "Rarefy10000")
#' methods <- c("unique_seq_alpha", "shannon_alpha", "gini_simpson_alpha")
#'
#' # Construct mock nested list for alpha diversity
#' raw.rarefy.res <- list()
#' for (r in reps) {
#'   raw.rarefy.res[[r]] <- list()
#'   for (s in samples) {
#'     raw.rarefy.res[[r]][[s]] <- list()
#'     for (d in depths) {
#'       raw.rarefy.res[[r]][[s]][[d]] <- list()
#'       for (m in methods) {
#'         raw.rarefy.res[[r]][[s]][[d]][[m]] <- round(runif(1, 0, 5), 2)
#'       }
#'     }
#'   }
#' }
#'
#' # Compute averaged matrices for each metric
#' result_list <- get_alpha_metric_matrix(raw.rarefy.res, metrics = methods)
#'
#' # View results for each alpha diversity metric
#' result_list$unique_seq_alpha
#' result_list$shannon_alpha
#' result_list$gini_simpson_alpha
#' @export
get_alpha_metric_matrix = function(x, metrics) {
  reps = names(x)
  samples = unique(unlist(lapply(reps, function(rep) names(x[[rep]]))))

  all_depths = unique(unlist(lapply(reps, function(rep) {
    lapply(x[[rep]], names)
  })))

  result_list = list()

  for (metric in metrics) {
    mat = matrix(NA, nrow = length(samples), ncol = length(all_depths),
                 dimnames = list(samples, all_depths))

    for (sample in samples) {
      for (depth in all_depths) {
        vals = sapply(reps, function(rep) {
          val = tryCatch(
            x[[rep]][[sample]][[depth]][[metric]],
            error = function(e) NA
          )
          if (is.null(val)) NA else val
        })
        mat[sample, depth] = if (all(is.na(vals))) NA else mean(vals, na.rm = TRUE)
      }
    }

    result_list[[metric]] = as.data.frame(mat)
  }

  return(result_list)
}
