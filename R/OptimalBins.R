#' Exact partial R^2 for L given Z via t-statistic (QR-based)
#'
#' Computes the exact partial coefficient of determination for the library-size
#' regressor \eqn{L} in the linear model \eqn{y ~ 1 + L + Z}, defined as
#' \deqn{R^2_{L \mid Z} = t_L^2 / (t_L^2 + \mathrm{df}),}
#' where \eqn{t_L} is the usual t-statistic for the coefficient of \eqn{L} and
#' \eqn{\mathrm{df}} is the residual degrees of freedom. The implementation uses
#' QR decomposition with column pivoting for numerical stability, and extracts
#' \eqn{\mathrm{Var}(\hat\beta_L)} from \eqn{R^{-1}R^{-T}} without forming
#' \eqn{(X'X)^{-1}} explicitly.
#'
#' @param y Numeric vector of length \eqn{n}: response.
#' @param L Numeric vector of length \eqn{n}: focal regressor whose partial
#'   \eqn{R^2} is reported (e.g., \eqn{\log_{10} N}). Must correspond row-wise to \code{y}.
#' @param Z Optional numeric matrix or data frame with \eqn{n} rows and \eqn{q}
#'   columns of adjustment covariates. If \code{NULL}, the model is \eqn{y ~ 1 + L}.
#'
#' @return A single \code{numeric} in \eqn{[0,1]}: the exact partial \eqn{R^2}
#'   for \eqn{L \mid Z}. Returns \code{NA_real_} if the design is singular, the
#'   residual degrees of freedom are nonpositive, or a variance component is not
#'   finite (the caller can penalize such bins as needed).
#'
#' @details
#' Let \eqn{X = [1, L, Z]} with column pivoting \eqn{X P = Q R}. The OLS
#' estimator is obtained via QR; the residual variance is \eqn{s^2}. The variance
#' of \eqn{\hat\beta_L} is taken from the diagonal of \eqn{(X'X)^{-1}}
#' computed as \eqn{R^{-1}R^{-T}} in the \emph{permuted} column space, then
#' mapped back to the original ordering using the pivot. The t-statistic
#' \eqn{t_L = \hat\beta_L / \mathrm{se}(\hat\beta_L)} yields
#' \eqn{R^2_{L \mid Z} = t_L^2 / (t_L^2 + \mathrm{df})}, which equals the
#' squared partial correlation when there is a single focal regressor \eqn{L}.
#'
#' @section Assumptions and input hygiene:
#' Inputs must be finite numerics and aligned by row. This function does not
#' remove \code{NA}; pre-filter or impute upstream. If \eqn{L} or any column of
#' \eqn{Z} is constant, or multicollinearity makes the design rank-deficient,
#' the function returns \code{NA_real_}.
#'
#' @section Complexity:
#' One QR solve per call; computational cost is \eqn{O(n p^2)} with
#' \eqn{p = 2 + q}.
#'
#' @import changepointGA
#' @importFrom GAReg gareg_knots
#' @examples
#' set.seed(1)
#' n <- 200
#' Z1 <- rnorm(n); Z2 <- rbinom(n, 1, 0.4)
#' L  <- rlnorm(n, meanlog = 8, sdlog = 0.5)
#' y  <- 3 + 0.2*log10(L) + 0.5*Z1 - 0.8*Z2 + rnorm(n, sd = 1)
#'
#' # Partial R^2 via this function
#' R2_fun <- .partial_R2_via_t(y, log10(L), cbind(Z1, Z2))
#' R2_fun
#'
#' # Cross-check against lm(): t^2 / (t^2 + df)
#' fit <- lm(y ~ log10(L) + Z1 + Z2)
#' summ <- summary(fit)
#' tL   <- summ$coefficients["log10(L)", "t value"]
#' df   <- fit$df.residual
#' R2_lm <- tL^2 / (tL^2 + df)
#' all.equal(R2_fun, as.numeric(R2_lm), tolerance = 1e-10)
#'
#' @keywords internal
#' @noRd
.partial_R2_via_t <- function(y, L, Z = NULL) {

  X <- if (is.null(Z)) cbind(1, L) else cbind(1, L, Z)

  qx <- qr(X)
  r  <- qx$rank
  p  <- ncol(X)
  if (r < p) return(NA_real_)

  beta <- qr.coef(qx, y)

  # Residuals and residual degrees of freedom
  res <- y - X %*% beta
  df  <- length(y) - r
  if (df <= 0) return(NA_real_)
  s2 <- sum(res^2) / df
  R      <- qr.R(qx)
  Rinvt  <- backsolve(R, diag(p))
  Vdiag_perm <- colSums(Rinvt^2)
  pos_L <- which(qx$pivot == 2L)
  if (length(pos_L) != 1L) return(NA_real_)

  # Standard error
  se_L <- sqrt(s2 * Vdiag_perm[pos_L])
  if (!is.finite(beta[2]) || !is.finite(se_L) || se_L <= 0) return(NA_real_)

  # t-statistic
  tval <- as.numeric(beta[2] / se_L)
  # exact partial R^2
  tval^2 / (tval^2 + df)
}

.mbrarefy_selected_indices <- function(fit = NULL, best_knots = NULL, Lgrid) {
  Lgrid <- as.numeric(Lgrid)
  if (!is.null(best_knots)) {
    best_knots <- sort(unique(as.numeric(best_knots)))
    idx <- match(best_knots, Lgrid)
    if (any(is.na(idx))) {
      idx <- vapply(best_knots, function(z) which.min(abs(Lgrid - z)), integer(1L))
    }
    idx <- sort(unique(idx[idx >= 2L & idx <= (length(Lgrid) - 1L)]))
    return(idx)
  }

  if (is.null(fit)) stop("Provide either 'fit' or 'best_knots'.", call. = FALSE)
  idx <- as.integer(fit@bestsol)
  sort(unique(idx[is.finite(idx) & idx >= 2L & idx <= (length(Lgrid) - 1L)]))
}

.mbrarefy_bin_definition <- function(Lgrid, idx = integer(0L)) {
  Lgrid <- as.numeric(Lgrid)
  idx <- sort(unique(as.integer(idx)))
  idx <- idx[idx >= 2L & idx <= (length(Lgrid) - 1L)]
  anchor_idx <- c(1L, idx)
  depths_op <- Lgrid[anchor_idx]
  BinCuts <- c(depths_op, Inf)
  list(anchor_idx = anchor_idx, depths_op = depths_op, BinCuts = BinCuts)
}

.mbrarefy_bin_id <- function(Lorig, BinCuts) {
  cut(
    as.numeric(Lorig),
    breaks = as.numeric(BinCuts),
    right = FALSE,
    include.lowest = TRUE,
    labels = FALSE
  )
}

.mbrarefy_bin_widths <- function(Lgrid, anchor_idx) {
  upper_idx <- c(anchor_idx[-1L], length(Lgrid))
  pmax(as.numeric(Lgrid[upper_idx]) - as.numeric(Lgrid[anchor_idx]),
       .Machine$double.eps)
}


#' Fixed-m GA objective: lower-bound–anchored partial R^2 across bins
#'
#' Objective function for GA-based knot selection with a fixed number of interior
#' knots. Given a rarefaction depth grid \code{x_unique} (length \eqn{M}),
#' precomputed alpha-diversity matrix \code{Y} (\eqn{n \times M}), and original
#' library sizes \code{Lorig} (length \eqn{n}), this function:
#' \itemize{
#'   \item decodes a chromosome carrying \eqn{m} interior knot \emph{indices}
#'         on the grid (sentinel \eqn{M+1} marks the end),
#'   \item forms \eqn{m+1} bins on \code{Lorig} using the grid cutpoints,
#'   \item anchors each bin at its \strong{lower bound} rarefaction depth,
#'   \item computes the exact binwise partial \eqn{R^2} of \eqn{L} (optionally
#'         adjusted for \code{Z}) via \code{.partial_R2_via_t},
#'   \item returns a weighted average of those binwise \eqn{R^2} values.
#' }
#' Smaller values indicate weaker within-bin association between alpha diversity
#' (at the bin’s anchor depth) and the original library size.
#'
#' @param knot_bin Numeric vector chromosome of the form
#'   \code{c(m, tau_1, ..., tau_m, sentinel = M+1, ...)}. Only the first \code{m}
#'   \emph{interior} indices \code{tau_k} appearing before the sentinel are used.
#' @param plen Ignored (kept for API symmetry with GAreg).
#' @param y Ignored by this objective (kept for API symmetry with GAreg).
#' @param x Numeric vector; if \code{x_unique} is missing, the grid is
#'   \code{sort(unique(x))}.
#' @param x_unique Strictly increasing numeric grid of rarefaction depths
#'   (length \eqn{M}). Columns of \code{Y} must correspond to this grid.
#' @param x_base Ignored (API symmetry).
#' @param fixedknots Integer scalar \eqn{m}: the number of interior knots to
#'   extract from \code{knot_bin}.
#' @param Y Numeric matrix \eqn{n \times M}: per-sample alpha diversity evaluated
#'   at each grid depth. Entry \code{Y[i,j]} should be \code{NA} when
#'   \code{x_unique[j] > Lorig[i]}.
#' @param Lorig Numeric vector (length \eqn{n}): original library sizes used for
#'   binning and as the predictor in the partial-\eqn{R^2} calculation.
#' @param Z Optional numeric matrix \eqn{n \times q} of covariates to partial out.
#'   If \code{NULL}, the model is \eqn{y ~ 1 + L}.
#' @param weight Bin weighting scheme in the aggregate: \code{"count"} (bin size),
#'   \code{"equal"} (uniform), or \code{"width"} (grid width on \code{x_unique}).
#' @param min_subjects Integer; minimum usable subjects per bin. Bins that do not
#'   meet this or fail numerically are penalized with \eqn{R^2=1}.
#' @param minDist Optional integer; minimum spacing between consecutive interior
#'   knot indices (in grid steps). If provided, chromosomes violating it are infeasible.
#' @param use_logL Logical; if \code{TRUE} (default) use \eqn{\log_{10}(L)} as the
#'   regressor in the partial-\eqn{R^2}; otherwise use raw \code{Lorig}.
#'
#' @return A single numeric: the weighted mean of binwise partial \eqn{R^2} values.
#'   Returns \code{Inf} for infeasible chromosomes (e.g., bad indices/spacing, size
#'   mismatches).
#'
#' @details
#' \strong{Chromosome decoding.} Let \eqn{M = } \code{length(x_unique)}. The function
#' reads the tail of \code{knot_bin} up to the first occurrence of \eqn{M+1} (sentinel),
#' keeps the first \eqn{m=\code{fixedknots}} interior indices in \eqn{\{2,\dots,M-1\}},
#' enforces uniqueness and optional \code{minDist}, and sorts them.
#'
#' \strong{Bins and anchors.} With interior indices \eqn{\tau_1<\cdots<\tau_m},
#' define operational lower-bound depths \code{x_unique[c(1, tau)]} and bin
#' boundaries \code{c(x_unique[c(1, tau)], Inf)}. Subjects are assigned to
#' half-open bins \eqn{[c_b, c_{b+1})}. The anchor column for bin \eqn{b} is
#' the bin's lower-bound grid index.
#'
#' \strong{Binwise metric.} For each bin, the response is \code{Y[, anchor]}, the
#' predictor is \code{log10(Lorig)} if \code{use_logL=TRUE} else \code{Lorig}, optionally
#' partialling out \code{Z}. The exact partial \eqn{R^2} of \eqn{L \mid Z} is computed
#' via \code{.partial_R2_via_t}. If the design is singular, \code{df <= 0}, or
#' \eqn{\mathrm{se}(\hat\beta_L)} is not finite, the bin contributes \eqn{R^2=1}.
#'
#' \strong{Aggregation.} The objective is the weighted mean of the per-bin \eqn{R^2}
#' with weights chosen by \code{weight}. This value is minimized by the GA.
#'
#' @seealso \code{.partial_R2_via_t} for the exact partial-\eqn{R^2} calculation.
#'
#' @examples
#' set.seed(1)
#' ## toy grid and data
#' x_unique <- seq(1000, 5000, by = 1000)  # M = 5
#' n <- 150
#' Lorig <- sample(1500:5000, n, replace = TRUE)
#' ids <- seq_len(n)
#' # build Y: alpha increases with depth and library size (for illustration)
#' Y <- outer(ids, x_unique, function(i,j) pmin(400, 0.02*Lorig[i] + 0.05*j)) +
#'      matrix(rnorm(n*length(x_unique), sd=3), n)
#' Y[ t(matrix(rep(x_unique, each=n), nrow=n)) > Lorig ] <- NA  # beyond-N -> NA
#'
#' # chromosome with m=2 interior knots at indices 3 and 4 (i.e., x=3000, 4000)
#' M <- length(x_unique)
#' knot_bin <- c(2, 3, 4, M+1)
#'
#' # requires .partial_R2_via_t() to be defined
#' fixBinRegObj(knot_bin,
#'              x = x_unique,
#'              x_unique = x_unique,
#'              fixedknots = 2L,
#'              Y = Y,
#'              Lorig = Lorig,
#'              Z = NULL,
#'              weight = "count",
#'              min_subjects = 15L,
#'              minDist = 1L,
#'              use_logL = TRUE)
#'
#' @export
fixBinRegObj <- function(
    knot_bin,
    plen = 0,
    y,
    x,
    x_unique,
    x_base = NULL,       # unused (API symmetry from GAreg)
    fixedknots,
    Y,
    Lorig,
    Z = NULL,
    weight = c("count","equal","width"),
    min_subjects = 8L,
    minDist = NULL,
    use_logL = TRUE
) {

  weight <- match.arg(weight)

  # grid & data checks
  if (missing(x_unique) || is.null(x_unique)) x_unique <- sort(unique(x))
  x_unique <- sort(unique(x_unique))
  M <- length(x_unique)
  if (!is.matrix(Y)) stop("Y must be an n x length(x_unique) matrix.")
  n <- NROW(Y)
  if (length(Lorig) != n || NCOL(Y) != M) return(Inf)

  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (NROW(Z) != n) return(Inf)
    q <- NCOL(Z)
  } else {
    q <- 0L
  }

  # transform L if requested (recommended)
  L_use_raw <- as.numeric(Lorig)
  L_use <- if (use_logL) log10(pmax(L_use_raw, 1)) else L_use_raw

  # ----- decode fixed-m interior indices (sentinel = M+1) -----
  m       <- as.integer(fixedknots)
  tail    <- as.integer(knot_bin[-1L])
  end_pos <- match(M + 1L, tail)
  if (is.na(end_pos)) return(Inf)

  cand <- tail[seq_len(end_pos - 1L)]
  cand <- cand[cand != 0L]
  if (length(cand) < m) return(Inf)
  idx <- sort(cand[seq_len(m)])

  # interiority & spacing
  Lb <- 2L; Ub <- M - 1L
  if (any(!is.finite(idx)) || any(idx < Lb) || any(idx > Ub)) return(Inf)
  if (anyDuplicated(idx)) return(Inf)
  if (!is.null(minDist) && length(idx) > 1L && any(diff(idx) < as.integer(minDist))) return(Inf)

  # ----- bins on the grid; lower-bound column = anchor_idx[b] -----
  bin_def    <- .mbrarefy_bin_definition(x_unique, idx)
  anchor_idx <- bin_def$anchor_idx
  B          <- length(anchor_idx)

  # map subjects to bins by original library size: [lower, next lower)
  bin_id <- .mbrarefy_bin_id(L_use_raw, bin_def$BinCuts)

  # bin widths for optional weighting
  widths <- .mbrarefy_bin_widths(x_unique, anchor_idx)

  # per-bin metric & counts (initialize with penalty 1)
  R2b <- rep(1.0, B)
  nb  <- integer(B)

  for (b in seq_len(B)) {
    rows_b <- which(bin_id == b)
    if (!length(rows_b)) next

    jLB <- anchor_idx[b]
    y_b <- Y[rows_b, jLB, drop = TRUE]
    L_b <- L_use[rows_b]

    ok <- is.finite(y_b) & is.finite(L_b)
    if (!is.null(Z)) {
      Zb <- Z[rows_b, , drop = FALSE]
      ok <- ok & apply(Zb, 1L, function(.) all(is.finite(.)))
    }

    y_b <- y_b[ok]; L_b <- L_b[ok]
    if (!is.null(Z)) Zb <- Zb[ok, , drop = FALSE] else Zb <- NULL
    nb[b] <- length(y_b)

    # need enough subjects for intercept + L (+ q covariates)
    if (nb[b] < max(min_subjects, (if (is.null(Z)) 2L else (2L + q)) + 1L)) next

    R2 <- .partial_R2_via_t(y_b, L_b, Zb)
    if (is.finite(R2)) R2b[b] <- R2 else R2b[b] <- 1.0
  }

  # weights
  wb <- switch(weight,
               count = nb,
               equal = rep(1, B),
               width = pmax(widths, .Machine$double.eps))
  if (sum(wb) <= 0) return(Inf)
  wb <- wb / sum(wb)

  # weighted aggregate (smaller is better)
  sum(wb * R2b)
}

#' Varying-m GA objective: lower-bound–anchored partial R^2 across bins
#'
#' Objective for GA-based, \emph{variable}-number-of-knots selection. Given a
#' rarefaction depth grid \code{x_unique} of length \eqn{M}, a precomputed
#' per-sample alpha-diversity matrix \code{Y} (\eqn{n \times M}), and the
#' original library sizes \code{Lorig} (length \eqn{n}), this function:
#' \itemize{
#'   \item decodes a chromosome carrying an arbitrary number of interior knot
#'         \emph{indices} (the first sentinel \eqn{M+1} marks the end),
#'   \item forms bins on \code{Lorig} using those grid cutpoints,
#'   \item anchors each bin at its \strong{lower bound} rarefaction depth,
#'   \item computes the exact binwise partial \eqn{R^2} of \eqn{L} (optionally
#'         adjusted for \code{Z}) via \code{.partial_R2_via_t},
#'   \item returns a weighted average of the binwise \eqn{R^2} values.
#' }
#' Smaller values indicate weaker within-bin association between alpha diversity
#' (evaluated at each bin’s lower-bound anchor) and the original library size.
#'
#' @param knot_bin Numeric vector chromosome. The function infers the set of
#'   interior knot indices by reading \emph{all} entries after the first element
#'   up to (but not including) the first occurrence of the sentinel \eqn{M+1}.
#'   Only interior grid indices in \code{2:(M-1)} are kept; duplicates are removed;
#'   remaining indices are sorted.
#' @param plen Ignored (kept for API symmetry with GA wrappers).
#' @param y Ignored by this objective (API symmetry with GA wrappers).
#' @param x Numeric vector. If \code{x_unique} is missing, the grid is
#'   \code{sort(unique(x))}.
#' @param x_unique Strictly increasing numeric grid of rarefaction depths
#'   (length \eqn{M}). Columns of \code{Y} must correspond to this grid.
#' @param x_base Ignored (API symmetry).
#' @param Y Numeric matrix \eqn{n \times M}: alpha diversity per sample at each
#'   grid depth. Convention: set \code{Y[i,j] <- NA} when \code{x_unique[j] > Lorig[i]}.
#' @param Lorig Numeric vector (length \eqn{n}): original library sizes used for
#'   binning and as the predictor in the partial-\eqn{R^2} calculation.
#' @param Z Optional numeric matrix \eqn{n \times q} of covariates to partial out.
#'   If \code{NULL}, the model is \eqn{y ~ 1 + L}.
#' @param weight Bin weights in the aggregate: \code{"count"} (bin size),
#'   \code{"equal"} (uniform), or \code{"width"} (grid width on \code{x_unique}).
#' @param min_subjects Integer; minimum usable subjects per bin. Bins that do not
#'   meet this or fail numerically are penalized with \eqn{R^2=1}.
#' @param minDist Optional integer; minimum spacing between consecutive interior
#'   knot indices (in \emph{grid steps}). If provided and violated, the chromosome
#'   is infeasible.
#' @param min_knots Integer; minimum number of selected interior knots. Use
#'   \code{min_knots = 1L} to exclude the one-bin solution in varying-\eqn{K}
#'   sensitivity analyses.
#' @param use_logL Logical; if \code{TRUE} (default) uses \eqn{\log_{10}(Lorig)}
#'   as the regressor in the partial-\eqn{R^2}; otherwise uses raw \code{Lorig}.
#'
#' @return A single numeric: the weighted mean of binwise partial \eqn{R^2}
#'   values. Returns \code{Inf} for infeasible chromosomes (e.g., missing sentinel,
#'   no valid interior indices, spacing violation, misaligned dimensions, or
#'   zero-sum weights).
#'
#' @details
#' \strong{Chromosome decoding (variable m).} Let \eqn{M = } \code{length(x_unique)}.
#' Read \code{knot_bin[-1]} up to the first \eqn{M+1} sentinel; keep unique interior
#' indices in \code{2:(M-1)}, sort them, and optionally enforce a minimum spacing
#' \code{minDist} in grid units. The resulting set defines \eqn{m} interior knots.
#'
#' \strong{Bins and anchors.} With interior indices \eqn{\tau_1<\cdots<\tau_m},
#' define operational lower-bound depths \code{x_unique[c(1, tau)]} and bin
#' boundaries \code{c(x_unique[c(1, tau)], Inf)}. Subjects are assigned to
#' half-open bins \eqn{[c_b, c_{b+1})}. The anchor column for bin \eqn{b} is
#' the bin's lower-bound grid index.
#'
#' \strong{Per-bin metric.} In each bin, the response is \code{Y[, anchor]}, the
#' predictor is \code{log10(Lorig)} if \code{use_logL=TRUE} else \code{Lorig}, and
#' covariates \code{Z} (if provided) are partialled out. The bin’s statistic is the
#' exact partial \eqn{R^2} of \eqn{L \mid Z} computed by \code{.partial_R2_via_t}.
#' If the design is singular, \code{df <= 0}, or \eqn{\mathrm{se}(\hat\beta_L)} is
#' not finite, the bin contributes \eqn{R^2=1}.
#'
#' \strong{Aggregation.} The objective minimized by the GA is the weighted mean
#' of per-bin \eqn{R^2} with weights chosen via \code{weight}. No explicit penalty
#' on the number of knots is included here; use \code{minDist} or GA controls to
#' regularize if needed.
#'
#' @seealso \code{.partial_R2_via_t} for the exact partial-\eqn{R^2} calculation.
#'
#' @examples
#' set.seed(42)
#' # grid and toy data
#' x_unique <- seq(1000, 5000, by = 1000)  # M = 5
#' n <- 120
#' Lorig <- sample(1200:5200, n, replace = TRUE)
#' Y <- outer(seq_len(n), x_unique, function(i, d) 0.015*Lorig[i] + 0.04*d) +
#'      matrix(rnorm(n*length(x_unique), sd = 2), n)
#' Y[ t(matrix(rep(x_unique, each=n), nrow=n)) > Lorig ] <- NA  # beyond-N -> NA
#'
#' # chromosome with variable m: interior indices {3}; sentinel M+1 ends the encoding
#' M <- length(x_unique)
#' knot_bin <- c(NA_real_, 3, M+1)  # first element ignored by this objective
#'
#' # requires .partial_R2_via_t() to be defined
#' varBinRegObj(knot_bin,
#'              x = x_unique, x_unique = x_unique,
#'              Y = Y, Lorig = Lorig,
#'              Z = NULL,
#'              weight = "count",
#'              min_subjects = 15L,
#'              minDist = 1L,
#'              use_logL = TRUE)
#'
#' @export
varBinRegObj <- function(
    knot_bin,
    plen = 0,
    y,
    x,
    x_unique,
    x_base = NULL,
    Y,
    Lorig,
    Z = NULL,
    weight = c("count","equal","width"),
    min_subjects = 8L,
    minDist = NULL,
    min_knots = 0L,
    use_logL = TRUE
) {
  weight <- match.arg(weight)

  # grid & data checks
  if (missing(x_unique) || is.null(x_unique)) x_unique <- sort(unique(x))
  x_unique <- sort(unique(x_unique))
  M <- length(x_unique)
  if (!is.matrix(Y)) stop("Y must be an n x M matrix.")
  n <- NROW(Y); if (length(Lorig) != n || NCOL(Y) != M) return(Inf)
  if (!is.null(Z)) { Z <- as.matrix(Z); if (NROW(Z) != n) return(Inf); q <- NCOL(Z) } else q <- 0L

  tail    <- as.integer(knot_bin[-1L])
  end_pos <- match(M + 1L, tail)
  if (is.na(end_pos)) return(Inf)

  idx <- sort(unique(as.integer(tail[seq_len(end_pos - 1L)])))
  idx <- idx[idx >= 2L & idx <= (M - 1L)]
  if (length(idx) < as.integer(min_knots)) return(Inf)
  # optional spacing (grid units)
  if (!is.null(minDist) && length(idx) > 1L && any(diff(idx) < as.integer(minDist))) return(Inf)

  bin_def    <- .mbrarefy_bin_definition(x_unique, idx)
  anchor_idx <- bin_def$anchor_idx
  B          <- length(anchor_idx)

  # bin assignment uses RAW library sizes; bins are [lower, next lower)
  bin_id <- .mbrarefy_bin_id(Lorig, bin_def$BinCuts)

  # weights
  widths <- .mbrarefy_bin_widths(x_unique, anchor_idx)
  wb <- switch(weight,
               count = as.numeric(table(factor(bin_id, levels = seq_len(B)))),
               equal = rep(1, B),
               width = pmax(widths, .Machine$double.eps))
  if (sum(wb) <= 0) return(Inf)
  wb <- wb / sum(wb)

  # model predictor (association on log10 L is standard; set use_logL=FALSE to use raw L)
  L_model <- if (use_logL) log10(pmax(Lorig, 1)) else as.numeric(Lorig)

  # per-bin metric (penalize failures with 1)
  R2b <- rep(1.0, B)
  for (b in seq_len(B)) {
    rows_b <- which(bin_id == b)
    if (!length(rows_b)) next

    jLB <- anchor_idx[b]                          # lower-bound anchor column
    y_b <- Y[rows_b, jLB, drop = TRUE]
    L_b <- L_model[rows_b]

    ok  <- is.finite(y_b) & is.finite(L_b)
    if (!is.null(Z)) { Zb <- Z[rows_b, , drop = FALSE]; ok <- ok & apply(Zb, 1L, function(.) all(is.finite(.))) }
    y_b <- y_b[ok]; L_b <- L_b[ok]; if (!is.null(Z)) Zb <- Zb[ok, , drop = FALSE] else Zb <- NULL

    # enforce minimal subjects for intercept + L (+ q covariates)
    if (length(y_b) < max(min_subjects, 2L + q + 1L)) { R2b[b] <- 1.0; next }

    R2 <- .partial_R2_via_t(y_b, L_b, Zb)
    R2b[b] <- if (is.finite(R2)) R2 else 1.0
  }

  # objective = weighted mean of binwise partial R^2 (smaller is better)
  sum(wb * R2b)
}

#' Select MBRarefy Library-Size Bins
#'
#' User-facing wrapper for data-adaptive MBRarefy cutpoint selection. The
#' function calls \code{GAReg::gareg_knots()} with either the fixed-\eqn{K}
#' objective \code{\link{fixBinRegObj}} or the varying-\eqn{K} objective
#' \code{\link{varBinRegObj}}, then returns selected cutpoints, operational
#' lower-bound depths, bin boundaries, assignments, and the fitted GA object.
#'
#' @param Lorig Numeric vector of original library sizes, one per sample.
#' @param Lgrid Numeric increasing vector of candidate rarefaction depths.
#' @param Y Numeric matrix of alpha-diversity values with samples in rows and
#'   depths in columns. \code{nrow(Y)} must equal \code{length(Lorig)} and
#'   \code{ncol(Y)} must equal \code{length(Lgrid)}.
#' @param mode Character; \code{"fixed"} for user-specified bin count or
#'   \code{"varying"} for data-adaptive bin count.
#' @param K Integer number of bins used when \code{mode = "fixed"}. The number
#'   of interior cutpoints is \code{K - 1}. Defaults to six bins.
#' @param Z Optional covariate matrix adjusted in the residual library-size
#'   objective.
#' @param weight Bin weighting scheme passed to the objective:
#'   \code{"count"}, \code{"equal"}, or \code{"width"}.
#' @param min_subjects Minimum usable subjects per bin in the objective.
#' @param minDist Optional minimum spacing, in grid index units, between
#'   selected interior cutpoints.
#' @param use_logL Logical; if \code{TRUE}, the objective uses
#'   \eqn{\log_{10}(Lorig)} as the library-size regressor.
#' @param min_knots Minimum number of interior cutpoints for
#'   \code{mode = "varying"}. The final simulation setting uses
#'   \code{min_knots = 1L}.
#' @param gaMethod GA method passed to \code{GAReg::gareg_knots()}.
#' @param cptgactrl Optional control object from \code{GAReg::cptgaControl()}.
#'   If \code{NULL}, a higher-budget default is used.
#' @param seed Optional random seed for reproducibility.
#' @param ... Additional arguments passed to \code{GAReg::gareg_knots()}.
#'
#' @return A list with selected cutpoints, bin boundaries, assignments,
#'   operational lower-bound depths, the fitted GA object, and settings.
#'
#' @examples
#' \dontrun{
#' fit_bins <- selectMBRarefyBins(
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   Y = as.matrix(USC),
#'   mode = "fixed",
#'   K = 6L,
#'   min_subjects = 20L
#' )
#'
#' y_anchor <- extractMBRarefyAlpha(
#'   Y = as.matrix(USC),
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   BinCuts = fit_bins$BinCuts,
#'   depths_op = fit_bins$depths_op
#' )
#' }
#' @export
selectMBRarefyBins <- function(
    Lorig,
    Lgrid,
    Y,
    mode = c("fixed", "varying"),
    K = 6L,
    Z = NULL,
    weight = c("count", "equal", "width"),
    min_subjects = 20L,
    minDist = 1L,
    use_logL = TRUE,
    min_knots = 1L,
    gaMethod = "cptga",
    cptgactrl = NULL,
    seed = NULL,
    ...
) {
  mode <- match.arg(mode)
  weight <- match.arg(weight)
  Lorig <- as.numeric(Lorig)
  Lgrid <- as.numeric(Lgrid)
  Y <- as.matrix(Y)

  if (length(Lorig) != nrow(Y)) stop("length(Lorig) must equal nrow(Y).", call. = FALSE)
  if (length(Lgrid) != ncol(Y)) stop("length(Lgrid) must equal ncol(Y).", call. = FALSE)
  if (length(Lgrid) < 3L) stop("Lgrid must contain at least three depths.", call. = FALSE)
  if (any(!is.finite(Lgrid)) || any(diff(Lgrid) <= 0)) {
    stop("Lgrid must be finite and strictly increasing.", call. = FALSE)
  }
  if (any(!is.finite(Lorig))) stop("Lorig must contain finite values.", call. = FALSE)

  if (is.null(cptgactrl)) {
    cptgactrl <- GAReg::cptgaControl(
      popSize = 400,
      maxgen = 100000,
      pchangepoint = 0.3
    )
  }

  if (!is.null(seed)) set.seed(seed)

  common_args <- list(
    y = rep(0, length(Lgrid)),
    x = Lgrid,
    gaMethod = gaMethod,
    cptgactrl = cptgactrl,
    minDist = minDist,
    Y = Y,
    Lorig = Lorig,
    Z = Z,
    weight = weight,
    min_subjects = min_subjects,
    use_logL = use_logL,
    ...
  )

  if (mode == "fixed") {
    K <- as.integer(K)
    if (length(K) != 1L || !is.finite(K) || K < 2L) {
      stop("K must be a single integer at least 2 for fixed mode.", call. = FALSE)
    }
    if ((K - 1L) > (length(Lgrid) - 2L)) {
      stop("K is too large for the number of available interior grid points.", call. = FALSE)
    }
    fit <- do.call(
      GAReg::gareg_knots,
      c(common_args, list(ObjFunc = fixBinRegObj, fixedknots = K - 1L))
    )
  } else {
    fit <- do.call(
      GAReg::gareg_knots,
      c(common_args, list(ObjFunc = varBinRegObj, min_knots = min_knots))
    )
  }

  idx <- .mbrarefy_selected_indices(fit = fit, Lgrid = Lgrid)
  bin_def <- .mbrarefy_bin_definition(Lgrid, idx)
  bin_id <- .mbrarefy_bin_id(Lorig, bin_def$BinCuts)

  list(
    fit = fit,
    mode = mode,
    best_idx = idx,
    best_knots = Lgrid[idx],
    depths_op = bin_def$depths_op,
    BinCuts = bin_def$BinCuts,
    bin_id = bin_id,
    K = length(bin_def$depths_op),
    settings = list(
      weight = weight,
      min_subjects = min_subjects,
      minDist = minDist,
      use_logL = use_logL,
      min_knots = if (mode == "varying") min_knots else NA_integer_
    )
  )
}

#' Extract Bin-Anchored Alpha Diversity Values
#'
#' Extract the alpha-diversity value used for MBRarefy downstream inference:
#' each sample is assigned to a library-size bin and evaluated at that bin's
#' lower-bound rarefaction depth.
#'
#' @param Y Numeric sample-by-depth alpha-diversity matrix.
#' @param Lorig Numeric vector of original library sizes.
#' @param Lgrid Numeric vector of candidate rarefaction depths.
#' @param fit Optional fitted result returned by \code{selectMBRarefyBins()} or
#'   \code{GAReg::gareg_knots()}.
#' @param best_knots Optional numeric cutpoints used when \code{fit} is absent.
#' @param BinCuts Optional bin boundaries, usually from
#'   \code{selectMBRarefyBins()}.
#' @param depths_op Optional operational lower-bound depths, usually from
#'   \code{selectMBRarefyBins()}.
#' @param return Character; \code{"vector"} returns only anchored alpha values,
#'   while \code{"data.frame"} also returns bin IDs and rarefaction depths.
#'
#' @return A numeric vector or a data frame, depending on \code{return}.
#' @export
extractMBRarefyAlpha <- function(
    Y,
    Lorig,
    Lgrid,
    fit = NULL,
    best_knots = NULL,
    BinCuts = NULL,
    depths_op = NULL,
    return = c("vector", "data.frame")
) {
  return <- match.arg(return)
  Y <- as.matrix(Y)
  Lorig <- as.numeric(Lorig)
  Lgrid <- as.numeric(Lgrid)

  if (length(Lorig) != nrow(Y)) stop("length(Lorig) must equal nrow(Y).", call. = FALSE)
  if (length(Lgrid) != ncol(Y)) stop("length(Lgrid) must equal ncol(Y).", call. = FALSE)
  if (any(!is.finite(Lgrid)) || any(diff(Lgrid) <= 0)) {
    stop("Lgrid must be finite and strictly increasing.", call. = FALSE)
  }

  if (is.null(BinCuts) || is.null(depths_op)) {
    if (is.list(fit) && !is.null(fit$BinCuts) && !is.null(fit$depths_op)) {
      BinCuts <- fit$BinCuts
      depths_op <- fit$depths_op
    } else {
      idx <- .mbrarefy_selected_indices(fit = fit, best_knots = best_knots, Lgrid = Lgrid)
      bin_def <- .mbrarefy_bin_definition(Lgrid, idx)
      BinCuts <- bin_def$BinCuts
      depths_op <- bin_def$depths_op
    }
  }

  bin_id <- .mbrarefy_bin_id(Lorig, BinCuts)
  rarefy_depth <- depths_op[bin_id]
  col_id <- match(rarefy_depth, Lgrid)

  y_anchor <- rep(NA_real_, length(Lorig))
  ok <- !is.na(bin_id) & !is.na(col_id)
  y_anchor[ok] <- Y[cbind(which(ok), col_id[ok])]

  if (return == "vector") return(y_anchor)

  data.frame(
    alpha = y_anchor,
    bin_id = bin_id,
    rarefy_depth = rarefy_depth,
    library_size = Lorig
  )
}
