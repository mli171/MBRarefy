.mbrarefy_partial_R2_nested <- function(y, L, Z = NULL) {

  y <- as.numeric(y)
  L <- as.numeric(L)

  if (!is.null(Z)) Z <- as.matrix(Z)

  ok <- is.finite(y) & is.finite(L)

  if (!is.null(Z)) {ok <- ok & apply(Z, 1L, function(z) all(is.finite(z)))}

  y <- y[ok]
  L <- L[ok]

  if (!is.null(Z)) {Z <- Z[ok, , drop = FALSE]}

  n <- length(y)

  if (n < 4L) {return(NA_real_)}

  X0 <- if (is.null(Z)) {
    matrix(1, nrow = n, ncol = 1L)
  } else {
    cbind(1, Z)
  }

  X1 <- cbind(X0, L)

  q0 <- qr(X0)
  q1 <- qr(X1)

  if (q0$rank < ncol(X0) || q1$rank < ncol(X1)) {return(NA_real_)}

  SSE0 <- sum((y - qr.fitted(q0, y))^2)
  SSE1 <- sum((y - qr.fitted(q1, y))^2)

  if (!is.finite(SSE0) || !is.finite(SSE1) || SSE0 <= 0 || SSE1 < 0) {
    return(NA_real_)
  }

  SSE1 <- min(SSE1, SSE0)
  R2 <- 1 - SSE1 / SSE0

  min(max(R2, 0), 1)
}

.mbrarefy_fragmentation_penalty <- function(bin_id, B, min_subjects, n_total, eta_F = 0.05) {

  nb <- tabulate(bin_id, nbins = B)

  if (any(nb < min_subjects)) {return(Inf)}

  n_F <- eta_F * n_total
  F_raw <- mean(n_F / nb)
  F0 <- n_F / n_total
  F <- (F_raw - F0) / (1 - F0)

  min(max(F, 0), 1)
}

.mbrarefy_bin_definition <- function(Lgrid, idx = integer(0L)) {

  Lgrid <- as.numeric(Lgrid)
  idx <- sort(unique(as.integer(idx)))

  idx <- idx[
    idx >= 2L &
      idx <= (length(Lgrid) - 1L)
  ]

  anchor_idx <- c(1L, idx)
  depths_op <- Lgrid[anchor_idx]
  BinCuts <- c(depths_op, Inf)

  list(
    anchor_idx = anchor_idx,
    depths_op = depths_op,
    BinCuts = BinCuts
  )
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

.mbrarefy_evaluate_design <- function(idx, Lgrid, Lorig, Y, Z = NULL, min_subjects = 20L, eta_F = 0.05, return_details = FALSE) {

  Lgrid <- as.numeric(Lgrid)
  Lorig <- as.numeric(Lorig)
  Y <- as.matrix(Y)

  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (nrow(Z) != length(Lorig)) {
      stop("Z must have one row per observation in Lorig.", call. = FALSE)
    }
  }

  bin_def <- .mbrarefy_bin_definition(Lgrid, idx)
  bin_id <- .mbrarefy_bin_id(Lorig, bin_def$BinCuts)

  B <- length(bin_def$anchor_idx)
  nb <- tabulate(bin_id, nbins = B)

  ## Hard constraint
  if (any(nb < min_subjects)) {
    if (return_details) {
      return(list(valid = FALSE, objective = Inf))
    }
    return(Inf)
  }

  Lmod <- log10(pmax(Lorig, 1))
  Db <- numeric(B)

  for (b in seq_len(B)) {

    rows <- which(bin_id == b)
    j <- bin_def$anchor_idx[b]

    Zb <- if (is.null(Z)) {
      NULL
    } else {
      Z[rows, , drop = FALSE]
    }

    Db[b] <- .mbrarefy_partial_R2_nested(y = Y[rows, j], L = Lmod[rows], Z = Zb)

    if (!is.finite(Db[b])) {
      if (return_details) {return(list(valid = FALSE, objective = Inf))}
      return(Inf)
    }
  }

  ## Sample-size-weighted residual depth dependence
  Q <- sum((nb / sum(nb)) * Db)
  F <- .mbrarefy_fragmentation_penalty(bin_id = bin_id, B = B, min_subjects = min_subjects, n_total = length(Lorig), eta_F = eta_F)
  J <- Q + F

  if (!return_details) {return(J)}

  list(valid = TRUE, objective = J, Q = Q, F = F, K = length(idx),
       knot_index = idx, knot_depth = Lgrid[idx], bin_n = nb,
       lower_depth = bin_def$depths_op, bin_id = bin_id)
}




#' Objective Function for Fixed-K MBRarefy Cutpoint Selection
#'
#' Evaluates a candidate fixed-K library-size partition for
#' data-adaptive multi-bin rarefying.
#'
#' The number of interior cutpoints, \eqn{K}, is fixed in advance, and
#' the optimizer searches over their locations on the candidate
#' rarefaction-depth grid.
#'
#' Candidate partitions are evaluated using
#'
#' \deqn{J(\tau) = Q(\tau) + F(\tau),}
#'
#' where \eqn{Q} is the sample-size-weighted average of bin-specific
#' partial \eqn{R^2} values measuring residual library-size dependence,
#' and \eqn{F} is the sample-fragmentation penalty.
#'
#' Partitions containing invalid cutpoint configurations, cutpoints
#' outside the interior of the candidate grid, cutpoints violating the
#' minimum-spacing requirement, or bins with fewer than
#' \code{min_subjects} samples are assigned an infinite objective value.
#'
#' This function is primarily intended for use as the objective function
#' in the fixed-K optimizer used by \code{selectMBRarefyBins()}.
#'
#' @param knot_bin Candidate chromosome supplied by
#'   \code{GAReg::gareg_knots()}.
#' @param plen Compatibility argument used by the genetic-algorithm
#'   optimizer. It is not used directly by this objective function.
#' @param y Compatibility argument supplied by the optimizer. It is not
#'   used directly by this objective function.
#' @param x Numeric vector of candidate rarefaction depths.
#' @param x_unique Optional vector of unique candidate rarefaction
#'   depths. If missing, it is constructed from \code{x}.
#' @param x_base Compatibility argument supplied by the optimizer. It is
#'   not used directly by this objective function.
#' @param fixedknots Positive integer giving the fixed number of interior
#'   cutpoints.
#' @param Y Numeric matrix of alpha-diversity values with samples in rows
#'   and candidate rarefaction depths in columns.
#' @param Lorig Numeric vector of original library sizes.
#' @param Z Optional matrix of nuisance covariates used when quantifying
#'   residual library-size dependence.
#' @param min_subjects Positive integer giving the minimum number of
#'   samples required in every selected bin.
#' @param minDist Optional minimum spacing between interior cutpoints in
#'   grid-index units.
#' @param eta_F Scaling parameter for the fragmentation penalty.
#'
#' @return A numeric objective value. Smaller values indicate a preferred
#'   partition. Invalid or infeasible partitions return \code{Inf}.
#'
#' @seealso \code{\link{selectMBRarefyBins}}, \code{\link{varBinRegObj}}
#' @export
fixBinRegObj <- function(knot_bin, plen = 0, y, x, x_unique, x_base = NULL, fixedknots, Y, Lorig, Z = NULL, min_subjects = 20L, minDist = NULL, eta_F = 0.05) {

  if (missing(x_unique) || is.null(x_unique)) {
    x_unique <- sort(unique(x))
  }

  x_unique <- sort(unique(x_unique))
  M <- length(x_unique)

  m <- as.integer(fixedknots)

  tail <- as.integer(knot_bin[-1L])
  end_pos <- match(M + 1L, tail)

  if (is.na(end_pos)) return(Inf)

  cand <- tail[seq_len(end_pos - 1L)]
  cand <- cand[cand != 0L]

  if (length(cand) < m) return(Inf)

  idx <- sort(cand[seq_len(m)])

  if (length(idx) != m || anyDuplicated(idx) || any(idx < 2L) || any(idx > M - 1L)) {
    return(Inf)
  }

  if (!is.null(minDist) && length(idx) > 1L && any(diff(idx) < as.integer(minDist))) {
    return(Inf)
  }

  .mbrarefy_evaluate_design(
    idx = idx,
    Lgrid = x_unique,
    Lorig = Lorig,
    Y = Y,
    Z = Z,
    min_subjects = min_subjects,
    eta_F = eta_F,
    return_details = FALSE
  )
}





.mbrarefy_parse_joint_chromosome <- function(chromosome, M) {

  chromosome <- as.integer(chromosome)

  if (length(chromosome) < 2L) { return(NULL)}

  K <- chromosome[1L]

  if (!is.finite(K) || K < 0L || length(chromosome) < K + 2L) {return(NULL)}

  knots <- if (K == 0L) {
    integer(0)
  } else {
    sort(chromosome[2L:(K + 1L)])
  }

  if (length(knots) != K || anyDuplicated(knots) || (K > 0L && (any(knots < 2L) || any(knots > M - 1L)))) {
    return(NULL)
  }

  list(K = K, knots = knots)
}






#' Objective Function for Varying-K MBRarefy Cutpoint Selection
#'
#' Evaluates a candidate varying-K library-size partition for
#' data-adaptive multi-bin rarefying.
#'
#' The chromosome jointly encodes the number and locations of interior
#' cutpoints. Its first entry gives the number of interior cutpoints,
#' \eqn{K}, followed by the corresponding cutpoint indices.
#'
#' Candidate partitions are evaluated using
#'
#' \deqn{J(\tau) = Q(\tau) + F(\tau),}
#'
#' where \eqn{Q} is the sample-size-weighted average of bin-specific
#' partial \eqn{R^2} values measuring residual library-size dependence
#' and \eqn{F} is the sample-fragmentation penalty.
#'
#' Partitions containing bins with fewer than \code{min_subjects}
#' samples, invalid cutpoint configurations, or numbers of cutpoints
#' outside the allowed range are assigned an infinite objective value.
#'
#' This function is primarily intended for use as the objective function
#' in the varying-K optimizer used by \code{selectMBRarefyBins()}.
#'
#' @param chromosome Integer chromosome encoding the candidate partition.
#'   The first entry is the number of interior cutpoints \eqn{K}; the
#'   following \eqn{K} entries are their grid indices.
#' @param plen Compatibility argument used by the genetic-algorithm
#'   optimizer. It is not used directly by this objective function.
#' @param depth_grid Strictly increasing numeric vector of candidate
#'   rarefaction depths.
#' @param Lorig Numeric vector of original library sizes.
#' @param Y Numeric matrix of alpha-diversity values with samples in rows
#'   and candidate rarefaction depths in columns.
#' @param Z Optional matrix of nuisance covariates used when quantifying
#'   residual library-size dependence.
#' @param min_subjects Positive integer giving the minimum number of
#'   samples required in every selected bin.
#' @param eta_F Scaling parameter for the fragmentation penalty.
#' @param Kmin_allowed Minimum allowed number of interior cutpoints.
#'   The default is 1.
#' @param Kmax_allowed Maximum allowed number of interior cutpoints.
#'
#' @return A numeric objective value. Smaller values indicate a preferred
#'   partition. Invalid or infeasible partitions return \code{Inf}.
#'
#' @seealso \code{\link{selectMBRarefyBins}}, \code{\link{fixBinRegObj}}
#' @export
varBinRegObj <- function(chromosome, plen = 0, depth_grid, Lorig, Y, Z = NULL, min_subjects = 20L, eta_F = 0.05, Kmin_allowed = 1L, Kmax_allowed = Inf) {

  design <- .mbrarefy_parse_joint_chromosome(chromosome, length(depth_grid))

  if (is.null(design) || design$K < Kmin_allowed || design$K > Kmax_allowed) {
    return(Inf)
  }

  .mbrarefy_evaluate_design(
    idx = design$knots,
    Lgrid = depth_grid,
    Lorig = Lorig,
    Y = Y,
    Z = Z,
    min_subjects = min_subjects,
    eta_F = eta_F,
    return_details = FALSE
  )
}



.mbrarefy_make_ga_suggestions <- function(N_ga, Kmax) {

  out <- vector("list", Kmax + 1L)
  out[[1L]] <- NULL

  if (Kmax >= 1L) {
    for (K in seq_len(Kmax)) {
      pos <- round(seq(1, N_ga + 1, length.out = K + 2L))[2L:(K + 1L)]
      pos <- pmax(2L, pmin(N_ga, pos))
      out[[K + 1L]] <- if (length(unique(pos)) == K) {
        pos
      } else {
        NULL
      }
    }
  }

  out
}

#' Select MBRarefy Library-Size Cutpoints
#'
#' Select library-size cutpoints for data-adaptive multi-bin rarefying.
#'
#' The primary fixed-\eqn{K} mode prespecifies the number of interior
#' cutpoints and optimizes their locations. The optional varying-\eqn{K}
#' mode jointly selects the number and locations of interior cutpoints.
#'
#' Candidate partitions are evaluated using
#'
#' \deqn{J(\tau) = Q(\tau) + F(\tau),}
#'
#' where \eqn{Q} is the sample-size-weighted average of bin-specific
#' partial \eqn{R^2} values measuring residual library-size dependence,
#' and \eqn{F} is a fragmentation penalty that discourages unnecessary
#' or highly imbalanced partitioning.
#'
#' The criterion targets the magnitude of residual library-size
#' dependence rather than the statistical significance of a separate
#' within-bin test. Thus, selected partitions are not required to make
#' the library-size association non-significant in every individual bin.
#' Additional cutpoints are favored only when their reduction in
#' residual depth dependence is sufficient to offset the accompanying
#' increase in sample fragmentation.
#'
#' Every selected bin must contain at least \code{min_subjects} samples.
#' In this function, \eqn{K} always denotes the number of interior
#' cutpoints, so the corresponding number of bins is \eqn{B = K + 1}.
#' Under varying-\eqn{K} selection, \eqn{K >= 1}; therefore, the
#' simplest allowed partition contains one interior cutpoint and two bins.
#'
#' @param Lorig Numeric vector of original library sizes.
#' @param Lgrid Strictly increasing numeric vector of candidate
#'   rarefaction depths.
#' @param Y Numeric matrix of alpha-diversity values with samples in rows
#'   and candidate rarefaction depths in columns.
#' @param mode Character string specifying \code{"fixed"} or
#'   \code{"varying"} cutpoint selection.
#' @param K Positive integer giving the number of interior cutpoints for
#'   fixed-\eqn{K} selection. The resulting partition contains
#'   \code{K + 1} bins.
#' @param Kmax Positive integer giving the maximum number of interior
#'   cutpoints considered under varying-\eqn{K} selection.
#' @param Z Optional matrix of nuisance covariates used when quantifying
#'   residual library-size dependence. The primary exposure or group
#'   variable used for downstream association analysis should generally
#'   not be included when the partition is intended to be selected
#'   independently of that association.
#' @param min_subjects Positive integer giving the minimum number of
#'   samples required in every selected bin.
#' @param eta_F Scaling parameter for the fragmentation penalty.
#' @param minDist Minimum spacing between interior cutpoints in
#'   grid-index units.
#' @param fixed_control Optional control object returned by
#'   \code{GAReg::cptgaControl()} for fixed-\eqn{K} optimization.
#' @param varying_popSize Population size per island for varying-\eqn{K}
#'   optimization.
#' @param varying_numIslands Number of islands.
#' @param varying_maxMig Maximum number of migrations.
#' @param varying_maxgen Number of generations per migration.
#' @param varying_maxconv Convergence limit.
#' @param varying_pcrossover Crossover probability.
#' @param varying_pmutation Mutation probability.
#' @param varying_pchangepoint Changepoint-mutation probability.
#' @param seed Optional random seed for reproducible optimization.
#'
#' @return A list containing the fitted optimizer, selected cutpoint
#'   indices and depths, bin definitions, objective value and components,
#'   bin sample sizes, and optimization settings. The returned
#'   \code{K} is the number of interior cutpoints and \code{B} is the
#'   number of bins.
#'
#' @examples
#' \dontrun{
#' ## Fixed-K selection: three interior cutpoints, giving four bins
#' fit_fixed <- selectMBRarefyBins(
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   Y = as.matrix(USC),
#'   mode = "fixed",
#'   K = 3L,
#'   min_subjects = 20L,
#'   eta_F = 0.05,
#'   seed = 1
#' )
#'
#' fit_fixed$K
#' fit_fixed$B
#' fit_fixed$best_knots
#' fit_fixed$Q
#' fit_fixed$F
#' fit_fixed$objective
#'
#' ## Varying-K selection: jointly choose the number and locations
#' ## of interior cutpoints
#' fit_varying <- selectMBRarefyBins(
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   Y = as.matrix(USC),
#'   mode = "varying",
#'   Kmax = 10L,
#'   min_subjects = 20L,
#'   eta_F = 0.05,
#'   seed = 1
#' )
#'
#' fit_varying$K
#' fit_varying$B
#' fit_varying$best_knots
#' fit_varying$Q
#' fit_varying$F
#' fit_varying$objective
#'
#' ## Extract bin-anchored alpha-diversity values
#' y_anchor <- extractMBRarefyAlpha(
#'   Y = as.matrix(USC),
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   fit = fit_varying
#' )
#' }
#' @export
selectMBRarefyBins <- function(Lorig, Lgrid, Y, mode = c("fixed", "varying"), K = 3L, Kmax = 10L, Z = NULL, min_subjects = 20L,
                               eta_F = 0.05, minDist = 1L, fixed_control = NULL,
                               varying_popSize = 100L, varying_numIslands = 5L, varying_maxMig = 200L, varying_maxgen = 100L,
                               varying_maxconv = 201L, varying_pcrossover = 0.95, varying_pmutation = 0.30, varying_pchangepoint = 0.10,
                               seed = NULL) {

  mode <- match.arg(mode)

  Lorig <- as.numeric(Lorig)
  Lgrid <- as.numeric(Lgrid)
  Y <- as.matrix(Y)

  if (length(Lorig) != nrow(Y)) {
    stop("length(Lorig) must equal nrow(Y).", call. = FALSE)
  }

  if (length(Lgrid) != ncol(Y)) {
    stop("length(Lgrid) must equal ncol(Y).", call. = FALSE)
  }

  if (length(Lgrid) < 3L) {
    stop("Lgrid must contain at least three candidate depths.", call. = FALSE)
  }

  if (any(!is.finite(Lgrid)) || any(diff(Lgrid) <= 0)) {
    stop("Lgrid must be finite and strictly increasing.", call. = FALSE)
  }

  if (any(!is.finite(Lorig))) {
    stop("Lorig must contain only finite values.", call. = FALSE)
  }

  if (any(Lorig < Lgrid[1L])) {
    stop("All library sizes must be at least the minimum depth in Lgrid.", call. = FALSE)
  }

  if (!is.null(Z)) {
    Z <- as.matrix(Z)

    if (nrow(Z) != length(Lorig)) {
      stop("Z must have one row per sample.", call. = FALSE)
    }
  }

  min_subjects <- as.integer(min_subjects)

  if (length(min_subjects) != 1L || !is.finite(min_subjects) || min_subjects < 1L) {
    stop("min_subjects must be a positive integer.", call. = FALSE)
  }

  if (length(eta_F) != 1L || !is.finite(eta_F) || eta_F <= 0 || eta_F >= 1) {
    stop("eta_F must lie in (0, 1).", call. = FALSE)
  }

  minDist <- as.integer(minDist)

  if (length(minDist) != 1L || !is.finite(minDist) || minDist < 1L) {
    stop("minDist must be a positive integer.", call. = FALSE)
  }

  if (!is.null(seed)) {set.seed(seed)}

  M <- length(Lgrid)
  n <- length(Lorig)

  # FIEXED-K #

  if (mode == "fixed") {

    K <- as.integer(K)

    if (length(K) != 1L || !is.finite(K) || K < 1L) {
      stop("K must be a positive integer giving the number of interior cutpoints.", call. = FALSE)
    }

    if (K > M - 2L) {
      stop("K exceeds the number of available interior grid locations.", call. = FALSE)
    }

    if (n < (K + 1L) * min_subjects) {
      stop("The requested K is infeasible under the current min_subjects constraint.", call. = FALSE)
    }

    if (is.null(fixed_control)) {
      fixed_control <- GAReg::cptgaControl(popSize = 100L, maxgen = 7500L, pchangepoint = 0.10)
    }

    fit <- GAReg::gareg_knots(
      y = rep(0, M),
      x = Lgrid,
      ObjFunc = fixBinRegObj,
      minDist = minDist,
      gaMethod = "cptga",
      cptgactrl = fixed_control,
      fixedknots = K,
      Y = Y,
      Lorig = Lorig,
      Z = Z,
      min_subjects = min_subjects,
      eta_F = eta_F
    )

    idx <- sort(as.integer(fit@bestsol))

    components <- .mbrarefy_evaluate_design(
      idx = idx,
      Lgrid = Lgrid,
      Lorig = Lorig,
      Y = Y,
      Z = Z,
      min_subjects = min_subjects,
      eta_F = eta_F,
      return_details = TRUE
    )

    if (!isTRUE(components$valid)) {
      stop("Fixed-K optimization did not return a valid partition.", call. = FALSE)
    }

    bin_def <- .mbrarefy_bin_definition(Lgrid, idx)

    return(
      list(
        fit = fit,
        mode = "fixed",
        K = K,
        B = K + 1L,
        best_idx = idx,
        best_knots = Lgrid[idx],
        depths_op = bin_def$depths_op,
        BinCuts = bin_def$BinCuts,
        bin_id = components$bin_id,
        objective = components$objective,
        Q = components$Q,
        F = components$F,
        bin_n = components$bin_n,
        components = components,
        settings = list(
          min_subjects = min_subjects,
          eta_F = eta_F,
          minDist = minDist,
          seed = seed
        )
      )
    )
  }


  # VARYING-K #

  Kmax <- as.integer(Kmax)

  if (length(Kmax) != 1L || !is.finite(Kmax) || Kmax < 1L) {
    stop("Kmax must be a positive integer.", call. = FALSE)
  }

  N_ga <- M - 1L

  Kmax_use <- min(Kmax, floor(n / min_subjects) - 1L, N_ga - 1L)

  if (Kmax_use < 1L) {
    stop("No feasible positive K is available under the current min_subjects constraint.", call. = FALSE)
  }

  lmax_safe <- N_ga + 2L

  fit <- changepointGA::cptgaisl(
    ObjFunc = varBinRegObj,
    N = N_ga,
    popSize = as.integer(varying_popSize),
    numIslands = as.integer(varying_numIslands),
    pcrossover = varying_pcrossover,
    pmutation = varying_pmutation,
    pchangepoint = varying_pchangepoint,
    minDist = minDist,
    mmax = Kmax_use,
    lmax = lmax_safe,
    maxMig = as.integer(varying_maxMig),
    maxgen = as.integer(varying_maxgen),
    maxconv = as.integer(varying_maxconv),
    option = "cp",
    monitoring = FALSE,
    parallel = FALSE,
    tol = 1e-8,
    seed = seed,
    suggestions = .mbrarefy_make_ga_suggestions(N_ga = N_ga, Kmax = Kmax_use),
    depth_grid = Lgrid,
    Lorig = Lorig,
    Y = Y,
    Z = Z,
    min_subjects = min_subjects,
    eta_F = eta_F,
    Kmin_allowed = 1L,
    Kmax_allowed = Kmax_use
  )

  chromosome <- fit@overbestchrom

  design <- .mbrarefy_parse_joint_chromosome(chromosome = chromosome, M = M)

  if (is.null(design)) {
    stop("Could not parse the best varying-K chromosome.", call. = FALSE)
  }

  components <- .mbrarefy_evaluate_design(
    idx = design$knots,
    Lgrid = Lgrid,
    Lorig = Lorig,
    Y = Y,
    Z = Z,
    min_subjects = min_subjects,
    eta_F = eta_F,
    return_details = TRUE
  )

  if (!isTRUE(components$valid)) {
    stop("Varying-K optimization did not return a valid partition.", call. = FALSE)
  }

  bin_def <- .mbrarefy_bin_definition(Lgrid, design$knots)

  list(
    fit = fit,
    mode = "varying",
    K = design$K,
    B = design$K + 1L,
    best_idx = design$knots,
    best_knots = Lgrid[design$knots],
    depths_op = bin_def$depths_op,
    BinCuts = bin_def$BinCuts,
    bin_id = components$bin_id,
    objective = components$objective,
    Q = components$Q,
    F = components$F,
    bin_n = components$bin_n,
    chromosome = chromosome,
    components = components,
    settings = list(
      Kmax = Kmax,
      Kmax_used = Kmax_use,
      min_subjects = min_subjects,
      eta_F = eta_F,
      minDist = minDist,
      seed = seed,
      popSize = varying_popSize,
      numIslands = varying_numIslands,
      maxMig = varying_maxMig,
      maxgen = varying_maxgen,
      maxconv = varying_maxconv,
      pcrossover = varying_pcrossover,
      pmutation = varying_pmutation,
      pchangepoint = varying_pchangepoint
    )
  )
}



#' Extract Bin-Anchored Alpha-Diversity Values
#'
#' Extract the alpha-diversity value used for downstream MBRarefy
#' association analysis. Each sample is assigned to a selected
#' library-size bin and evaluated at that bin's lower-bound
#' rarefaction depth.
#'
#' The partition can be supplied using a fitted object returned by
#' \code{selectMBRarefyBins()}, a vector of interior cutpoint depths,
#' or explicit bin boundaries and operational rarefaction depths.
#'
#' @param Y Numeric sample-by-depth alpha-diversity matrix, with samples
#'   in rows and candidate rarefaction depths in columns.
#' @param Lorig Numeric vector of original library sizes.
#' @param Lgrid Strictly increasing numeric vector of candidate
#'   rarefaction depths.
#' @param fit Optional result returned by \code{selectMBRarefyBins()}.
#' @param best_knots Optional numeric vector of interior cutpoint depths
#'   when \code{fit} is not supplied.
#' @param BinCuts Optional numeric vector of bin boundaries. When supplied
#'   directly, \code{depths_op} must also be supplied.
#' @param depths_op Optional numeric vector of bin-specific lower-bound
#'   rarefaction depths.
#' @param return Character string; \code{"vector"} returns only the
#'   bin-anchored alpha-diversity values, while \code{"data.frame"}
#'   additionally returns bin assignments, rarefaction depths, and
#'   original library sizes.
#'
#' @return A numeric vector or data frame, depending on \code{return}.
#'
#' @examples
#' \dontrun{
#' fit_bins <- selectMBRarefyBins(
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   Y = as.matrix(USC),
#'   mode = "fixed",
#'   K = 3L,
#'   seed = 1
#' )
#'
#' y_anchor <- extractMBRarefyAlpha(
#'   Y = as.matrix(USC),
#'   Lorig = dataPheno$totalReads,
#'   Lgrid = depths,
#'   fit = fit_bins
#' )
#' }
#'
#' @export
extractMBRarefyAlpha <- function(Y, Lorig, Lgrid, fit = NULL,
                                 best_knots = NULL, BinCuts = NULL, depths_op = NULL,
                                 return = c("vector", "data.frame")) {

  return <- match.arg(return)

  Y <- as.matrix(Y)
  Lorig <- as.numeric(Lorig)
  Lgrid <- as.numeric(Lgrid)

  if (length(Lorig) != nrow(Y)) {
    stop("length(Lorig) must equal nrow(Y).", call. = FALSE)
  }

  if (length(Lgrid) != ncol(Y)) {
    stop("length(Lgrid) must equal ncol(Y).", call. = FALSE)
  }

  if (any(!is.finite(Lgrid)) || any(diff(Lgrid) <= 0)) {
    stop("Lgrid must be finite and strictly increasing.", call. = FALSE)
  }

  if (any(!is.finite(Lorig))) {
    stop("Lorig must contain only finite values.", call. = FALSE)
  }

  if (any(Lorig < Lgrid[1L])) {
    stop("All library sizes must be at least the minimum depth in Lgrid.", call. = FALSE)
  }

  if (!is.null(fit)) {
    if (!is.list(fit) || is.null(fit$BinCuts) || is.null(fit$depths_op)) {
      stop("fit must be an object returned by selectMBRarefyBins().", call. = FALSE)
    }
    BinCuts <- fit$BinCuts
    depths_op <- fit$depths_op
  } else if (!is.null(best_knots)) {
    best_knots <- sort(as.numeric(best_knots))
    idx <- match(best_knots, Lgrid)
    if (anyNA(idx)) {
      stop("All best_knots must be present in Lgrid.", call. = FALSE)
    }
    bin_def <- .mbrarefy_bin_definition(Lgrid = Lgrid, idx = idx)
    BinCuts <- bin_def$BinCuts
    depths_op <- bin_def$depths_op
  } else {
    if (is.null(BinCuts) || is.null(depths_op)) {
      stop(paste0("Supply either fit, best_knots, or both ", "BinCuts and depths_op."), call. = FALSE)
    }
    BinCuts <- as.numeric(BinCuts)
    depths_op <- as.numeric(depths_op)
  }

  BinCuts <- as.numeric(BinCuts)
  depths_op <- as.numeric(depths_op)

  if (length(BinCuts) != length(depths_op) + 1L) {
    stop("length(BinCuts) must equal length(depths_op) + 1.", call. = FALSE)
  }

  if (!is.infinite(BinCuts[length(BinCuts)])) {
    stop("The final entry of BinCuts must be Inf.", call. = FALSE)
  }

  if (any(diff(BinCuts) <= 0)) {
    stop("BinCuts must be strictly increasing.", call. = FALSE)
  }

  col_id <- match(depths_op, Lgrid)

  if (anyNA(col_id)) {
    stop("All depths_op values must be present in Lgrid.", call. = FALSE)
  }

  bin_id <- .mbrarefy_bin_id(Lorig = Lorig, BinCuts = BinCuts)

  if (anyNA(bin_id)) {
    stop("Some samples could not be assigned to a library-size bin.", call. = FALSE)
  }

  rarefy_depth <- depths_op[bin_id]
  sample_col <- col_id[bin_id]
  y_anchor <- Y[cbind(seq_along(Lorig), sample_col)]

  if (return == "vector") {return(as.numeric(y_anchor))}

  data.frame(alpha = as.numeric(y_anchor), bin_id = bin_id, rarefy_depth = rarefy_depth, library_size = Lorig)

}



.mbrarefy_selected_indices <- function(fit = NULL, best_knots = NULL, Lgrid) {

  Lgrid <- as.numeric(Lgrid)

  ## Explicit cutpoint depths take priority
  if (!is.null(best_knots)) {
    best_knots <- sort(unique(as.numeric(best_knots)))
    idx <- match(best_knots, Lgrid)
    if (anyNA(idx)) {
      stop("All best_knots must be present in Lgrid.", call. = FALSE)
    }

    idx <- sort(unique(as.integer(idx)))

    if (any(idx < 2L) || any(idx > length(Lgrid) - 1L)) {
      stop("Selected cutpoints must be interior points of Lgrid.", call. = FALSE)
    }

    return(idx)
  }

  ## New MBRarefy wrapper object
  if (is.list(fit) && !is.null(fit$best_idx)) {

    idx <- sort(unique(as.integer(fit$best_idx)))

    if (any(idx < 2L) || any(idx > length(Lgrid) - 1L)) {
      stop("fit contains invalid cutpoint indices.", call. = FALSE)
    }

    return(idx)
  }

  ## Backward compatibility with GAReg fixed-K objects
  if (!is.null(fit) && isS4(fit) && "bestsol" %in% methods::slotNames(fit)) {
    idx <- sort(unique(as.integer(fit@bestsol)))
    idx <- idx[is.finite(idx) & idx >= 2L & idx <= length(Lgrid) - 1L]
    return(idx)
  }

  stop("Provide best_knots or a fitted object containing selected cutpoints.", call. = FALSE)

}
