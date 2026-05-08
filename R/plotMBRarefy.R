#' Plot MBRarefy-selected library-size bins
#'
#' Visualize library-size binning and alpha-diversity values used in
#' multibin rarefying analysis.
#'
#' This function produces a \code{ggplot2} scatterplot of alpha diversity
#' against original library size. It supports two plotting modes. In
#' \code{y_mode = "bin_lower_bound"}, each sample is plotted using the
#' alpha-diversity value evaluated at the lower-bound rarefaction depth of
#' its assigned MBRarefy bin. This corresponds to the bin-anchored diversity
#' used in MBRarefy downstream inference. In \code{y_mode = "fixed_depth"},
#' all samples are plotted using alpha diversity evaluated at a single
#' user-specified rarefaction depth, which is useful for illustrating
#' conventional single-depth rarefying.
#'
#' Selected MBRarefy cutpoints can be supplied either through a fitted
#' \code{GAReg::gareg_knots()} object via \code{fit}, or directly through
#' \code{best_knots}. When \code{fit} is supplied, the function extracts
#' the selected knot indices from \code{fit@bestsol} and maps them to the
#' corresponding depths in \code{Lgrid}.
#'
#' @param Lorig Numeric vector of original sample library sizes. Its length
#'   must equal the number of rows in \code{Y}.
#' @param Lgrid Numeric vector of candidate rarefaction depths. Its length
#'   must equal the number of columns in \code{Y}.
#' @param Y Numeric matrix or matrix-like object of alpha-diversity values,
#'   with samples in rows and rarefaction depths in columns. Entry
#'   \code{Y[i, j]} is the alpha diversity for sample \code{i} rarefied to
#'   depth \code{Lgrid[j]}. Missing values are allowed when a sample cannot
#'   be rarefied to a requested depth.
#' @param fit Optional fitted \code{GAReg::gareg_knots()} object containing
#'   selected cutpoint indices in \code{fit@bestsol}. Required when
#'   \code{best_knots} is not supplied and \code{y_mode = "bin_lower_bound"}.
#' @param best_knots Optional numeric vector of selected library-size
#'   cutpoints. If supplied, these are used directly instead of extracting
#'   cutpoints from \code{fit}.
#' @param y_mode Character string specifying the plotting mode. Choices are
#'   \code{"bin_lower_bound"} and \code{"fixed_depth"}. The
#'   \code{"bin_lower_bound"} mode plots each sample at its assigned
#'   bin-specific lower-bound rarefaction depth. The \code{"fixed_depth"}
#'   mode plots all samples at a single global rarefaction depth.
#' @param fixed_depth Numeric rarefaction depth used when
#'   \code{y_mode = "fixed_depth"}. The closest available value in
#'   \code{Lgrid} is used.
#' @param y_label Optional character string for the y-axis label. If
#'   \code{NULL}, a default label is generated based on \code{y_mode}.
#' @param title_prefix Optional character string for the plot title. If
#'   \code{NULL}, a default title is generated based on \code{y_mode}.
#' @param base_size Numeric base font size passed to
#'   \code{ggplot2::theme_bw()}.
#' @param numLegendRow Integer number of rows used for the color and shape
#'   legends.
#' @param show_missing_as_rug Logical; if \code{TRUE}, samples with
#'   unavailable alpha-diversity values at the plotted depth are shown as
#'   rug marks along the x-axis.
#' @param show_bin_coloring Logical; only used when
#'   \code{y_mode = "fixed_depth"}. If \code{TRUE} and bin cutpoints are
#'   available, points are colored by MBRarefy bin. If \code{FALSE}, all
#'   points are assigned to a single \code{"Samples"} group.
#' @param show_knot_lines Logical; if \code{TRUE}, vertical dashed blue
#'   lines are drawn at selected MBRarefy cutpoints when
#'   \code{y_mode = "bin_lower_bound"}.
#' @param show_fixed_depth_line Logical; if \code{TRUE}, a vertical dashed
#'   red line is drawn at \code{fixed_depth} when
#'   \code{y_mode = "fixed_depth"}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{plot}}{A \code{ggplot} object.}
#'   \item{\code{data}}{Data frame used to construct the plot.}
#'   \item{\code{data_valid}}{Subset of \code{data} with finite x and y values.}
#'   \item{\code{data_missing}}{Subset of \code{data} with finite library size but missing alpha-diversity values.}
#'   \item{\code{best_knots}}{Numeric vector of selected MBRarefy cutpoints.}
#'   \item{\code{depths_op}}{Operational rarefaction depths used as bin lower bounds.}
#'   \item{\code{BinCuts}}{Numeric vector of bin boundaries used for assigning samples.}
#' }
#'
#' @details
#' In \code{"bin_lower_bound"} mode, the function first defines operational
#' depths as \code{c(min(Lgrid), best_knots)} and assigns samples to bins
#' using \code{c(depths_op, Inf)} as bin boundaries. For each sample, the
#' plotted alpha-diversity value is extracted from the column of \code{Y}
#' corresponding to the lower-bound depth of that sample's assigned bin.
#' Samples with library sizes below the first operational depth are excluded
#' from this mode because their bin-anchored alpha diversity is undefined.
#'
#' In \code{"fixed_depth"} mode, the function selects the closest available
#' grid depth to \code{fixed_depth} and plots all samples using that single
#' rarefaction depth. This mode is useful for comparing MBRarefy with
#' conventional single-depth rarefying.
#'
#' @examples
#' \dontrun{
#' ## MBRarefy bin-lower-bound plot using a fitted fixed-K GA object
#' res_fix <- plotMBRarefy(
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   Y = as.matrix(USC),
#'   fit = fix_fit,
#'   y_mode = "bin_lower_bound",
#'   title_prefix = "Select cutpoints to reduce depth-diversity dependence (Fixed-K)",
#'   y_label = "Observed richness at bin lower bound",
#'   show_missing_as_rug = FALSE
#' )
#' res_fix$plot
#'
#' ## MBRarefy bin-lower-bound plot using a fitted varying-K GA object
#' res_var <- plotMBRarefy(
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   Y = as.matrix(USC),
#'   fit = var_fit,
#'   y_mode = "bin_lower_bound",
#'   title_prefix = "Select cutpoints to reduce depth-diversity dependence (Vary-K)",
#'   y_label = "Observed richness at bin lower bound",
#'   show_missing_as_rug = FALSE
#' )
#' res_var$plot
#'
#' ## Conventional single-depth rarefying diagnostic
#' res_40k <- plotMBRarefy(
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   Y = as.matrix(USC),
#'   y_mode = "fixed_depth",
#'   fixed_depth = 40000,
#'   show_bin_coloring = FALSE,
#'   show_knot_lines = FALSE,
#'   show_fixed_depth_line = TRUE
#' )
#' res_40k$plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_rug
#' @importFrom ggplot2 scale_shape_manual scale_x_log10 labs theme_bw theme
#' @importFrom ggplot2 element_text element_blank element_rect guide_legend guides
#' @importFrom ggplot2 expansion
#' @export
plotMBRarefy <- function(
    Lorig,
    Lgrid,
    Y,
    fit = NULL,
    best_knots = NULL,
    y_mode = c("bin_lower_bound", "fixed_depth"),
    fixed_depth = NULL,
    y_label = NULL,
    title_prefix = NULL,
    base_size = 11,
    numLegendRow = 1,
    show_missing_as_rug = TRUE,
    show_bin_coloring = TRUE,
    show_knot_lines = TRUE,
    show_fixed_depth_line = TRUE
) {

  y_mode <- match.arg(y_mode)

  Lorig <- as.numeric(Lorig)
  Lgrid <- as.numeric(Lgrid)
  Y <- as.matrix(Y)

  stopifnot(length(Lorig) == nrow(Y))
  stopifnot(length(Lgrid) == ncol(Y))
  stopifnot(all(is.finite(Lorig)))

  ## Extract selected cutpoints if available
  if (is.null(best_knots) && !is.null(fit)) {
    best_idx <- as.integer(fit@bestsol)
    best_idx <- sort(unique(best_idx[
      is.finite(best_idx) &
        best_idx >= 2L &
        best_idx <= (length(Lgrid) - 1L)
    ]))
    best_knots <- Lgrid[best_idx]
  }

  if (!is.null(best_knots)) {
    best_knots <- sort(unique(as.numeric(best_knots)))
    depths_op <- sort(unique(c(min(Lgrid), best_knots)))
    BinCuts <- c(depths_op, Inf)
  } else {
    depths_op <- min(Lgrid)
    BinCuts <- NULL
  }

  if (y_mode == "bin_lower_bound") {
    if (is.null(best_knots) && is.null(fit)) {
      stop("For y_mode = 'bin_lower_bound', provide 'fit' or 'best_knots'.")
    }

    keep <- is.finite(Lorig) & Lorig >= min(depths_op)

    L_use <- Lorig[keep]
    Y_use <- Y[keep, , drop = FALSE]

    bin_id <- cut(L_use, breaks = BinCuts, right = FALSE, include.lowest = TRUE, labels = FALSE)
    bin_lower_depth <- depths_op[bin_id]
    bins <- cut(L_use, breaks = BinCuts, right = FALSE, include.lowest = TRUE)

    col_id <- match(bin_lower_depth, Lgrid)
    y_plot <- Y_use[cbind(seq_len(nrow(Y_use)), col_id)]

    df_sc <- data.frame(N = L_use, Y = y_plot, Bin = bins, RarefyDepth = bin_lower_depth)

    if (is.null(y_label)) {
      y_label <- "Alpha diversity at the lower bound of each bin"
    }

    if (is.null(title_prefix)) {
      title_prefix <- "Select cutpoints to reduce depth-diversity dependence"
    }
  }

  if (y_mode == "fixed_depth") {
    if (is.null(fixed_depth)) {
      stop("Provide fixed_depth when y_mode = 'fixed_depth'.")
    }

    d <- Lgrid[which.min(abs(Lgrid - fixed_depth))]
    col_id <- match(d, Lgrid)

    if (show_bin_coloring && !is.null(BinCuts)) {
      bins <- cut(Lorig, breaks = BinCuts, right = FALSE, include.lowest = TRUE)
    } else {
      bins <- factor(rep("Samples", length(Lorig)))
    }

    df_sc <- data.frame(N = Lorig, Y = Y[, col_id], Bin = bins, RarefyDepth = d)

    if (is.null(y_label)) {
      y_label <- paste0("Alpha diversity at ", d, " library size")
    }
    if (is.null(title_prefix)) {
      title_prefix <- paste0("Alpha diversity at fixed rarefaction depth = ", d)
    }
  }

  df_valid <- df_sc[is.finite(df_sc$N) & is.finite(df_sc$Y), , drop = FALSE]
  df_missing <- df_sc[is.finite(df_sc$N) & !is.finite(df_sc$Y), , drop = FALSE]

  shape_pool <- c(16, 17, 15, 18, 3, 4, 8, 9, 0, 1, 2, 5, 6, 7)
  shape_vals <- shape_pool[seq_len(nlevels(df_valid$Bin))]

  x_log_breaks <- seq(
    floor(log10(min(df_valid$N, na.rm = TRUE)) * 10) / 10,
    ceiling(log10(max(df_valid$N, na.rm = TRUE)) * 10) / 10,
    by = 0.2
  )

  p <- ggplot2::ggplot(
    df_valid,
    ggplot2::aes(x = N, y = Y, shape = Bin, color = Bin)
  ) +
    ggplot2::geom_point(alpha = 0.6, stroke = 1) +
    ggplot2::scale_shape_manual(values = shape_vals) +
    ggplot2::scale_x_log10(
      breaks = 10^x_log_breaks,
      labels = x_log_breaks,
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(
      title = title_prefix,
      x = expression(log[10] * "(Library size)"),
      y = y_label
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1.2, fill = NA),
      axis.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = numLegendRow),
      shape = ggplot2::guide_legend(nrow = numLegendRow)
    )

  if (show_knot_lines && !is.null(best_knots) && y_mode == "bin_lower_bound") {
    p <- p + ggplot2::geom_vline(xintercept = best_knots, color = "blue", linetype = 2)
  }

  if (show_fixed_depth_line && y_mode == "fixed_depth") {
    p <- p + ggplot2::geom_vline(
      xintercept = fixed_depth,
      color = "red",
      linetype = 2
    )
  }

  if (show_missing_as_rug && nrow(df_missing) > 0) {
    p <- p +
      ggplot2::geom_rug(
        data = df_missing,
        ggplot2::aes(x = N),
        inherit.aes = FALSE,
        sides = "b",
        alpha = 0.7
      ) +
      ggplot2::labs(
        caption = paste0(
          nrow(df_missing),
          " samples have unavailable alpha diversity at the plotted depth and are shown as rug marks."
        )
      )
  }

  return(list(
    plot = p,
    data = df_sc,
    data_valid = df_valid,
    data_missing = df_missing,
    best_knots = best_knots,
    depths_op = depths_op,
    BinCuts = BinCuts
  ))
}


