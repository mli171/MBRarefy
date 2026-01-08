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
#' \strong{Bins and anchors.} With interior indices \eqn{\tau_1<\cdots<\tau_m}, define
#' \code{edges_idx <- c(1, tau, M)} and cutpoints \code{breaks <- x_unique[edges_idx]}.
#' Subjects are assigned by \eqn{(breaks_b, breaks_{b+1}]} on \code{Lorig}. The anchor
#' column for bin \eqn{b} is \code{edges_idx[b]} (the bin’s \emph{lower} grid bound).
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
#' @importFrom GAReg gareg_knots
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
  if (!is.null(minDist) && any(diff(idx) <= as.integer(minDist))) return(Inf)

  # ----- bins on the grid; lower-bound column = edges_idx[b] -----
  edges_idx <- c(1L, idx, M)
  breaks    <- x_unique[edges_idx]
  B         <- length(edges_idx) - 1L

  # map subjects to bins by original library size ( (breaks[b], breaks[b+1]] )
  bin_id <- findInterval(L_use_raw, breaks, rightmost.closed = TRUE, all.inside = TRUE)

  # bin widths for optional weighting
  widths <- x_unique[edges_idx[-1L]] - x_unique[edges_idx[-length(edges_idx)]]

  # per-bin metric & counts (initialize with penalty 1)
  R2b <- rep(1.0, B)
  nb  <- integer(B)

  for (b in seq_len(B)) {
    rows_b <- which(bin_id == b)
    if (!length(rows_b)) next

    jLB <- edges_idx[b]
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
#' set \code{edges_idx <- c(1, tau, M)} and \code{breaks <- x_unique[edges_idx]}.
#' Assign samples by \eqn{(breaks_b,\,breaks_{b+1}]} on \code{Lorig}. The anchor
#' column for bin \eqn{b} is \code{edges_idx[b]} (the bin’s \emph{lower} grid bound).
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
#' @importFrom GAReg gareg_knots
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
  # optional spacing (grid units)
  if (!is.null(minDist) && length(idx) > 1L && any(diff(idx) <= as.integer(minDist))) return(Inf)

  edges_idx <- c(1L, idx, M)
  breaks    <- x_unique[edges_idx]
  B         <- length(edges_idx) - 1L

  # bin assignment uses RAW library sizes; bins are (breaks[b], breaks[b+1]]
  bin_id <- findInterval(Lorig, breaks, rightmost.closed = TRUE, all.inside = TRUE)

  # weights
  widths <- x_unique[edges_idx[-1L]] - x_unique[edges_idx[-length(edges_idx)]]
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

    jLB <- edges_idx[b]                           # lower-bound anchor column
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
