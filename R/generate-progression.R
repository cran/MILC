# Tumor growth and disease progression (model_description/model.tex,
# Section 5): Gompertz growth from one malignant cell and three random volume
# thresholds (regional, distant, diagnosis), drawn exactly from their joint
# truncated log-normal law, by row-wise rejection (default, faster) or by
# conditional inverse-CDF draws (three uniforms per person; for common random
# numbers in calibration).

# Volume of one malignant cell and the maximum volume (mm^3), archival t_prog.
MILC_V0 <- pi * 0.01^3 / 6
MILC_VMAX <- pi * 130^3 / 6
# Bounds of log Z, with Z = V / 1000 the volume in ml.
MILC_LOGZ_MIN <- log(MILC_V0 / 1000)
MILC_LOGZ_MAX <- log(MILC_VMAX / 1000)

#' Time from onset to a tumor volume (Gompertz growth)
#'
#' Vectorized archival `t_vol`: \eqn{\tau(v) = -(365 m)^{-1}\log(1 -
#' \log(v/V_0)/31)} years, with \eqn{V_0 = \pi 0.01^3/6} mm\eqn{^3} and shape
#' \eqn{s = 31m}.
#'
#' @param vol volume (mm\eqn{^3}).
#' @param m Gompertz scale parameter (per day).
#' @return years from onset.
#' @keywords internal
milc_time_to_volume <- function(vol, m) {
  -log1p(-log(vol / MILC_V0) / 31) / (365 * m)
}

#' Draw volumes from a log-normal truncated to \eqn{(V_0, V_{max})}
#'
#' Volume \eqn{V = 1000 Z} mm\eqn{^3}, \eqn{Z \sim LN(\mu, \sigma)} in ml,
#' truncated to \eqn{V_0 < V < V_{max}}, by inverse CDF (exact).
#'
#' @param meanlog,sdlog log-normal parameters of \eqn{Z}, one per draw.
#' @return volumes (mm\eqn{^3}).
#' @keywords internal
draw_truncated_volume <- function(meanlog, sdlog) {
  n <- length(meanlog)
  lo <- stats::pnorm((log(MILC_V0 / 1000) - meanlog) / sdlog)
  hi <- stats::pnorm((log(MILC_VMAX / 1000) - meanlog) / sdlog)
  p <- lo + stats::runif(n) * (hi - lo)
  v <- 1000 * exp(meanlog + sdlog * stats::qnorm(p))
  # guard against rounding onto the bounds
  pmin(pmax(v, MILC_V0 * (1 + 1e-12)), MILC_VMAX * (1 - 1e-12))
}

#' Gauss-Legendre nodes and weights on \eqn{[-1, 1]}
#'
#' Golub-Welsch: eigen-decomposition of the Jacobi matrix.
#'
#' @param k number of nodes.
#' @return list with `x` (nodes, increasing) and `w` (weights).
#' @keywords internal
gauss_legendre <- function(k) {
  j <- seq_len(k - 1L)
  b <- j / sqrt(4 * j^2 - 1)
  J <- matrix(0, k, k)
  J[cbind(j, j + 1L)] <- b
  J[cbind(j + 1L, j)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
}

GL10 <- gauss_legendre(10L)

#' Normal probability of an interval, accurate in both tails
#'
#' \eqn{\Phi((hi - \mu)/\sigma) - \Phi((lo - \mu)/\sigma)}, computed from
#' upper-tail probabilities when the interval lies above the mean.
#'
#' @param lo,hi interval bounds (recycled).
#' @param mean,sd normal parameters.
#' @return probabilities, same shape as `hi` when `hi` is a matrix.
#' @keywords internal
pnorm_interval <- function(lo, hi, mean, sd) {
  zl <- (lo - mean) / sd
  zh <- (hi - mean) / sd
  out <- stats::pnorm(zh) - stats::pnorm(zl)
  upper <- zl > 0
  if (any(upper)) {
    zl <- rep_len(zl, length(zh))
    out[upper] <- stats::pnorm(zl[upper], lower.tail = FALSE) -
      stats::pnorm(zh[upper], lower.tail = FALSE)
  }
  out
}

#' Quantile function of a truncated normal
#'
#' Inverse CDF of \eqn{N(\mu, \sigma^2)} truncated to \eqn{(lo, hi)}, using
#' upper-tail probabilities when \eqn{lo > \mu} so that both tails are
#' accurate. Vectorized; arguments are recycled to the longest.
#'
#' @param u probabilities in \eqn{(0, 1)}.
#' @param mean,sd normal parameters.
#' @param lo,hi truncation bounds.
#' @return quantiles in \eqn{[lo, hi]}.
#' @keywords internal
qtruncnorm <- function(u, mean, sd, lo, hi) {
  n <- max(length(u), length(mean), length(sd), length(lo), length(hi))
  u <- rep_len(u, n)
  mean <- rep_len(mean, n)
  sd <- rep_len(sd, n)
  zl <- (rep_len(lo, n) - mean) / sd
  zh <- (rep_len(hi, n) - mean) / sd
  z <- numeric(n)
  upper <- zl > 0
  if (any(upper)) {
    i <- which(upper)
    pl <- stats::pnorm(zl[i], lower.tail = FALSE)
    ph <- stats::pnorm(zh[i], lower.tail = FALSE)
    z[i] <- stats::qnorm(pl - u[i] * (pl - ph), lower.tail = FALSE)
  }
  if (any(!upper)) {
    i <- which(!upper)
    pl <- stats::pnorm(zl[i])
    ph <- stats::pnorm(zh[i])
    z[i] <- stats::qnorm(pl + u[i] * (ph - pl))
  }
  mean + sd * pmin(pmax(z, zl), zh)
}

#' Marginal law of the larger of two ordered truncated normals
#'
#' For independent \eqn{X \sim N(\mu_X, \sigma_X^2)} and \eqn{Y \sim
#' N(\mu_Y, \sigma_Y^2)} restricted to \eqn{a < Y < X < b}, the marginal
#' density of \eqn{X} is proportional to \eqn{g(x) = \phi_X(x)\{\Phi_Y(x) -
#' \Phi_Y(a)\}} on \eqn{(a, b)}; its CDF is a bivariate normal probability with
#' no closed-form inverse. `ordered_marginal_table()` tabulates
#' \eqn{I(x) = \int_a^x g} on a grid by 10-point Gauss-Legendre quadrature per
#' cell. `ordered_marginal_q()` inverts \eqn{I(x)/I(b)}: it locates the cell
#' and solves \eqn{I(x_j) + \int_{x_j}^x g = u I(b)} within it by Newton steps
#' (derivative \eqn{g}) kept inside a shrinking bracket, with bisection when a
#' step leaves the bracket.
#'
#' @param mx,sx,my,sy normal parameters of \eqn{X} and \eqn{Y} (scalars).
#' @param a,b bounds.
#' @param n_cells grid cells.
#' @param u probabilities in \eqn{(0, 1)}.
#' @param tab table from `ordered_marginal_table()`.
#' @param tol absolute tolerance on \eqn{x}.
#' @param max_iter iteration cap.
#' @return `ordered_marginal_table()`: list with the parameters, `grid`, and
#'   `I` (cumulative integral at the grid points). `ordered_marginal_q()`:
#'   quantiles of \eqn{X}.
#' @keywords internal
ordered_marginal_table <- function(mx, sx, my, sy, a = MILC_LOGZ_MIN,
                                   b = MILC_LOGZ_MAX, n_cells = 4000L) {
  grid <- seq(a, b, length.out = n_cells + 1L)
  h <- diff(grid)
  mid <- (grid[-1L] + grid[-length(grid)]) / 2
  nodes <- outer(h / 2, GL10$x) + mid # cells x nodes
  g <- stats::dnorm(nodes, mx, sx) * pnorm_interval(a, nodes, my, sy)
  cell <- as.vector(g %*% GL10$w) * h / 2
  list(
    mx = mx, sx = sx, my = my, sy = sy, a = a, b = b,
    grid = grid, I = c(0, cumsum(cell))
  )
}

#' @rdname ordered_marginal_table
#' @keywords internal
ordered_marginal_q <- function(u, tab, tol = 1e-12, max_iter = 100L) {
  g_fun <- function(x) {
    stats::dnorm(x, tab$mx, tab$sx) * pnorm_interval(tab$a, x, tab$my, tab$sy)
  }
  # integral of g over [from, x], within one cell: 10-point Gauss-Legendre
  int_from <- function(from, x) {
    half <- (x - from) / 2
    nodes <- outer(half, GL10$x) + (from + x) / 2
    as.vector(g_fun(nodes) %*% GL10$w) * half
  }
  K <- length(tab$grid) - 1L
  target <- u * tab$I[K + 1L]
  j <- pmin(pmax(findInterval(target, tab$I), 1L), K)
  origin <- tab$grid[j] # integration origin, fixed
  rest <- target - tab$I[j] # mass to accumulate from the origin
  blo <- origin # bracket
  bhi <- tab$grid[j + 1L]
  w <- tab$I[j + 1L] - tab$I[j]
  x <- origin + (bhi - blo) * fifelse(w > 0, pmin(pmax(rest / w, 0), 1), 0.5)

  active <- seq_along(u)
  for (it in seq_len(max_iter)) {
    if (length(active) == 0L) break
    xa <- x[active]
    f <- int_from(origin[active], xa) - rest[active]
    above <- f > 0
    bhi[active[above]] <- xa[above]
    blo[active[!above]] <- xa[!above]
    x_new <- xa - f / g_fun(xa)
    out <- !is.finite(x_new) | x_new <= blo[active] | x_new >= bhi[active]
    x_new[out] <- (blo[active][out] + bhi[active][out]) / 2
    x[active] <- x_new
    done <- abs(x_new - xa) <= tol | f == 0 | (bhi[active] - blo[active]) <= tol
    active <- active[!done]
  }
  x
}

#' Volume thresholds by row-wise rejection
#'
#' \eqn{V_{diagn}} by inverse CDF of its truncated law; \eqn{V_{reg}} and
#' \eqn{V_{dist}} by inverse CDF of their laws truncated to
#' \eqn{(V_0, V_{max})}, redrawing only the rows with \eqn{V_{reg} \ge
#' V_{dist}}. Exact. Archival swap: \eqn{V_{reg}} uses `cdist`,
#' \eqn{V_{dist}} uses `creg`.
#'
#' @param p data.table of progression parameters, one row per draw.
#' @param max_rounds maximum rounds of redrawing.
#' @return list with `v_diagn`, `v_reg`, `v_dist` (mm\eqn{^3}).
#' @keywords internal
progression_volumes_rejection <- function(p, max_rounds = 1000L) {
  k <- nrow(p)
  v_diagn <- draw_truncated_volume(
    p$volume_diagnosis_meanlog, p$volume_diagnosis_sdlog
  )
  v_reg <- rep(NA_real_, k)
  v_dist <- rep(NA_real_, k)
  todo <- seq_len(k)
  for (round in seq_len(max_rounds)) {
    v_reg[todo] <- draw_truncated_volume(
      p$volume_distant_meanlog[todo], p$volume_distant_sdlog[todo]
    )
    v_dist[todo] <- draw_truncated_volume(
      p$volume_regional_meanlog[todo], p$volume_regional_sdlog[todo]
    )
    todo <- todo[v_reg[todo] >= v_dist[todo]]
    if (length(todo) == 0L) break
  }
  if (length(todo) > 0L) {
    stop("Order constraint V_reg < V_dist not met after ", max_rounds,
         " rounds for ", length(todo), " individual(s)")
  }
  list(v_diagn = v_diagn, v_reg = v_reg, v_dist = v_dist)
}

#' Volume thresholds by conditional inverse-CDF draws
#'
#' Three uniforms per person, no rejection: (i) \eqn{V_{diagn}} from its
#' log-normal truncated to \eqn{(V_0, V_{max})}; (ii) \eqn{V_{dist}} from its
#' marginal on the truncation set, density \eqn{\propto
#' f_{dist}(v)\{F_{reg}(v) - F_{reg}(V_0)\}}, by the numerically inverted CDF
#' of [ordered_marginal_table()]; (iii) \eqn{V_{reg}} given \eqn{V_{dist}} from
#' its log-normal truncated to \eqn{(V_0, V_{dist})}. Exact up to the
#' inversion tolerance (about \eqn{10^{-13}}). The map from the uniforms to the
#' volumes is deterministic and monotone in each uniform, so common random
#' numbers stay aligned across parameter values. About 5-6 times slower than
#' [progression_volumes_rejection()]. Archival swap as there.
#'
#' @inheritParams progression_volumes_rejection
#' @return list with `v_diagn`, `v_reg`, `v_dist` (mm\eqn{^3}).
#' @keywords internal
progression_volumes_conditional <- function(p) {
  k <- nrow(p)
  u_diagn <- stats::runif(k)
  u_dist <- stats::runif(k)
  u_reg <- stats::runif(k)
  lz_diagn <- qtruncnorm(
    u_diagn, p$volume_diagnosis_meanlog, p$volume_diagnosis_sdlog,
    MILC_LOGZ_MIN, MILC_LOGZ_MAX
  )
  lz_dist <- numeric(k)
  sets <- data.table(
    mx = p$volume_regional_meanlog, sx = p$volume_regional_sdlog,
    my = p$volume_distant_meanlog, sy = p$volume_distant_sdlog
  )
  sets[, row := .I]
  for (g in split(sets, by = c("mx", "sx", "my", "sy"))) {
    tab <- ordered_marginal_table(g$mx[1L], g$sx[1L], g$my[1L], g$sy[1L])
    lz_dist[g$row] <- ordered_marginal_q(u_dist[g$row], tab)
  }
  lz_reg <- qtruncnorm(
    u_reg, p$volume_distant_meanlog, p$volume_distant_sdlog,
    MILC_LOGZ_MIN, lz_dist
  )
  list(
    v_diagn = 1000 * exp(lz_diagn),
    v_reg = 1000 * exp(lz_reg),
    v_dist = 1000 * exp(lz_dist)
  )
}

#' Generate disease progression after onset
#'
#' Samples the volume thresholds \eqn{(V_{reg}, V_{dist}, V_{diagn})} from the
#' joint law of independent log-normals truncated to
#' \eqn{\{V_0 < V_{diagn} < V_{max}\} \cap \{V_0 < V_{reg} < V_{dist} <
#' V_{max}\}} (`model.tex`, Eq. truncset), exactly, by one of two methods:
#' `"rejection"` ([progression_volumes_rejection()]; default, faster) or
#' `"conditional"` ([progression_volumes_conditional()]; three uniforms per
#' person and a deterministic map from them to the volumes, for calibration
#' with common random numbers). The archival `t_prog` uses whole-vector
#' rejection.
#'
#' As in the archival `t_prog`, \eqn{V_{reg}} is drawn with the parameters
#' labelled `distant` (`cdist`) and \eqn{V_{dist}} with those labelled
#' `regional` (`creg`); see `model.tex`, Sections 5 and 9.
#'
#' @param age_onset age at onset of the first malignant cell; `NA` for no onset.
#' @param params data.frame with one row per person and columns
#'   `tumor_growth_rate` (\eqn{m}, per day), `volume_diagnosis_meanlog`,
#'   `volume_diagnosis_sdlog`, `volume_regional_meanlog`,
#'   `volume_regional_sdlog`, `volume_distant_meanlog`, `volume_distant_sdlog`
#'   (log-normal parameters of the volume in ml; archival `cdiagn`, `creg`,
#'   `cdist`).
#' @param method `"rejection"` or `"conditional"`; default
#'   `getOption("milc.progression_method")`.
#' @param max_rounds maximum rounds of redrawing (rejection only).
#' @return data.table with one row per person: `age_at_regional_disease`,
#'   `age_at_distal_disease`, `age_at_diagnosis`, `diameter_at_diagnosis` (mm),
#'   `stage_at_diagnosis` (`"localized"`, `"regional"`, `"distant"`); all `NA`
#'   where `age_onset` is `NA`. These are latent ages: diagnosis may fall after
#'   death.
#' @export
generate_progression <- function(age_onset, params,
                                 method = getOption("milc.progression_method"),
                                 max_rounds = 1000L) {
  method <- match.arg(method, c("rejection", "conditional"))
  n <- length(age_onset)
  stopifnot(nrow(params) == n)
  out <- data.table(
    age_at_regional_disease = rep(NA_real_, n),
    age_at_distal_disease = rep(NA_real_, n),
    age_at_diagnosis = rep(NA_real_, n),
    diameter_at_diagnosis = rep(NA_real_, n),
    stage_at_diagnosis = rep(NA_character_, n)
  )
  idx <- which(!is.na(age_onset))
  if (length(idx) == 0L) {
    return(out)
  }
  p <- as.data.table(params)[idx]
  vols <- switch(
    method,
    rejection = progression_volumes_rejection(p, max_rounds),
    conditional = progression_volumes_conditional(p)
  )
  v_diagn <- vols$v_diagn
  v_reg <- vols$v_reg
  v_dist <- vols$v_dist

  m <- p$tumor_growth_rate
  a <- age_onset[idx]
  set(out, idx, "age_at_regional_disease", a + milc_time_to_volume(v_reg, m))
  set(out, idx, "age_at_distal_disease", a + milc_time_to_volume(v_dist, m))
  set(out, idx, "age_at_diagnosis", a + milc_time_to_volume(v_diagn, m))
  set(out, idx, "diameter_at_diagnosis", (6 * v_diagn / pi)^(1 / 3))
  set(out, idx, "stage_at_diagnosis", fifelse(
    v_diagn < v_reg, "localized",
    fifelse(v_diagn < v_dist, "regional", "distant")
  ))
  out
}
