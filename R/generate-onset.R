# Onset of the first malignant cell: two-stage clonal expansion (TSCE) hazard,
# piecewise over the smoking history (model_description/model.tex, Section 4),
# sampled with nhppp::vdraw_cumulative_intensity() and a vectorized Newton
# inverse of the cumulative intensity (Section 2).

# Stable helpers ------------------------------------------------------------

# expm1(u) - u, without cancellation for small |u|.
expm1_minus_x <- function(u) {
  out <- expm1(u) - u
  small <- !is.na(u) & abs(u) < 0.1
  if (any(small)) {
    v <- u[small]
    out[small] <- v^2 * (1 / 2 + v * (1 / 6 + v * (1 / 24 + v * (1 / 120 +
      v * (1 / 720 + v * (1 / 5040 + v / 40320))))))
  }
  out
}

# log1p(x) - x, without cancellation for small |x|.
log1p_minus_x <- function(x) {
  out <- log1p(x) - x
  small <- !is.na(x) & abs(x) < 0.01
  if (any(small)) {
    v <- x[small]
    out[small] <- -v^2 * (1 / 2 - v * (1 / 3 - v * (1 / 4 - v * (1 / 5 -
      v * (1 / 6 - v / 7)))))
  }
  out
}

# TSCE kernel ----------------------------------------------------------------

#' TSCE constants for a parameter set
#'
#' For \eqn{\psi = (\nu, \alpha, \gamma, \mu)} and \eqn{X} stem cells
#' (`model.tex`, Section 4): \eqn{B = (-\gamma + \sqrt{\gamma^2 + 4\alpha\mu})/2},
#' computed as \eqn{2\alpha\mu / (\gamma + \sqrt{\gamma^2 + 4\alpha\mu})} to
#' avoid cancellation; \eqn{c = \gamma + 2B}; \eqn{K = \nu\mu X / (\gamma + B)};
#' and `nmx` \eqn{= \nu\mu X}. Vectorized.
#'
#' @param nu initiation rate of normal cells.
#' @param alpha division rate of initiated cells.
#' @param gamma \eqn{\alpha - \beta - \mu}.
#' @param mu malignant conversion rate.
#' @param X number of normal stem cells.
#' @return list with numeric vectors `B`, `c`, `K`, `nmx`.
#' @keywords internal
tsce_constants <- function(nu, alpha, gamma, mu, X) {
  root <- sqrt(gamma^2 + 4 * alpha * mu)
  B <- 2 * alpha * mu / (gamma + root)
  list(B = B, c = gamma + 2 * B, K = nu * mu * X / (gamma + B), nmx = nu * mu * X)
}

#' TSCE hazard and cumulative hazard
#'
#' `tsce_h()` is \eqn{h_\psi(t) = \nu\mu X (e^{ct}-1) / (c + B(e^{ct}-1))}.
#' `tsce_H()` is \eqn{\tilde H_\psi(t) = H_\psi(t) - H_\psi(0)}, where
#' \eqn{H_\psi} is the archival `HT_mal`, evaluated as
#' \eqn{K[(e^{ct}-1-ct)/c + \{\log(1+x) - x\}/B]} with \eqn{x = B(e^{ct}-1)/c},
#' which equals \eqn{K[-t + B^{-1}\log(1+x)]} without its cancellation.
#'
#' @param t ages (years), aligned with the elements of `k`.
#' @param k list from [tsce_constants()].
#' @return numeric vector.
#' @keywords internal
tsce_H <- function(t, k) {
  y <- expm1(k$c * t)
  k$K * (expm1_minus_x(k$c * t) / k$c + log1p_minus_x(k$B * y / k$c) / k$B)
}

#' @rdname tsce_H
#' @keywords internal
tsce_h <- function(t, k) {
  y <- expm1(k$c * t)
  k$nmx * y / (k$c + k$B * y)
}

subset_k <- function(k, idx) lapply(k, `[`, idx)

#' Per-person TSCE parameter sets
#'
#' Builds the non-smoking and smoking parameter sets
#' \eqn{\psi_{ns}} and \eqn{\psi_s(q)} of `model.tex`, Section 4, from the
#' columns of [deterministic_params] (archival `HT_mal` values):
#' \eqn{\kappa(q) = 1 + \alpha_1 q^{\alpha_2}}, \eqn{\alpha_s = \alpha_{s,0}
#' \kappa(q)}, \eqn{\gamma_s = \gamma\kappa(q)}, \eqn{\nu_s = \rho\nu},
#' \eqn{\mu_s = \mu}.
#'
#' @param params data.frame with one row per person and columns
#'   `n_stem_cells`, `normal_cell_initiation_rate`,
#'   `malignant_cell_conversion_rate`, `cell_division_rate_nonsmoking`,
#'   `cell_division_rate_smoking_base`, `net_proliferation_rate`,
#'   `smoking_initiation_multiplier`, `smoking_dose_coefficient`,
#'   `smoking_dose_exponent`.
#' @param intensity cigarettes per day while smoking (ignored for rows that
#'   never smoke; may be 0 or `NA` there).
#' @return list with elements `ns` and `s`, each a list from
#'   [tsce_constants()].
#' @keywords internal
tsce_person_constants <- function(params, intensity) {
  p <- params
  q <- intensity
  q[is.na(q)] <- 0
  kappa <- 1 + p$smoking_dose_coefficient * q^p$smoking_dose_exponent
  list(
    ns = tsce_constants(
      nu = p$normal_cell_initiation_rate,
      alpha = p$cell_division_rate_nonsmoking,
      gamma = p$net_proliferation_rate,
      mu = p$malignant_cell_conversion_rate,
      X = p$n_stem_cells
    ),
    s = tsce_constants(
      nu = p$smoking_initiation_multiplier * p$normal_cell_initiation_rate,
      alpha = p$cell_division_rate_smoking_base * kappa,
      gamma = p$net_proliferation_rate * kappa,
      mu = p$malignant_cell_conversion_rate,
      X = p$n_stem_cells
    )
  )
}

# Piecewise cumulative intensity and its inverse ---------------------------

#' Cumulative onset intensity over the smoking history
#'
#' \eqn{\Lambda_{mal}(t)} of `model.tex` Eq. (lambdaM): the non-smoking hazard
#' before `ts` and from `tq` on, the smoking hazard on `[ts, tq)`, all on the
#' attained-age clock. Never smokers have `ts = tq = Inf`; current smokers
#' `tq = Inf`.
#'
#' @param t numeric vector or matrix of ages; row `i` belongs to person `i`.
#' @param args list with `row_args`, a data.frame with columns `ts`, `tq`, and
#'   the per-person constants `B_ns, c_ns, K_ns, nmx_ns, B_s, c_s, K_s, nmx_s`
#'   (from [onset_row_args()]).
#' @return numeric vector or matrix like `t`.
#' @keywords internal
Lambda_mal <- function(t, args) {
  ra <- args$row_args
  if (is.matrix(t)) {
    out <- t
    for (j in seq_len(ncol(t))) out[, j] <- Lambda_mal(t[, j], args)
    return(out)
  }
  k_ns <- list(B = ra$B_ns, c = ra$c_ns, K = ra$K_ns, nmx = ra$nmx_ns)
  k_s <- list(B = ra$B_s, c = ra$c_s, K = ra$K_s, nmx = ra$nmx_s)
  ts <- ra$ts
  tq <- ra$tq
  t1 <- pmin(t, ts)
  out <- tsce_H(t1, k_ns)
  in_s <- !is.na(t) & t > ts
  if (any(in_s)) {
    i <- which(in_s)
    ks <- subset_k(k_s, i)
    out[i] <- out[i] + tsce_H(pmin(t[i], tq[i]), ks) - tsce_H(ts[i], ks)
  }
  after_q <- !is.na(t) & t > tq
  if (any(after_q)) {
    i <- which(after_q)
    kn <- subset_k(k_ns, i)
    out[i] <- out[i] + tsce_H(t[i], kn) - tsce_H(tq[i], kn)
  }
  out
}

#' Solve \eqn{\tilde H_\psi(t) = z} on \eqn{[lo, hi]} for each element
#'
#' Newton's method from the right end: \eqn{\tilde H_\psi} is increasing and
#' convex, so from a point with \eqn{\tilde H_\psi(t) \ge z} the iterates
#' decrease monotonically to the root. Elements not converged after `max_newton`
#' iterations (roots near \eqn{t = 0}, where \eqn{h \approx 0}) are finished by
#' bisection. Requires \eqn{\tilde H_\psi(lo) \le z \le \tilde H_\psi(hi)}.
#'
#' @param z targets.
#' @param lo,hi bracket.
#' @param k list from [tsce_constants()], aligned with `z`.
#' @param tol relative tolerance on `t`.
#' @param max_newton,max_bisect iteration caps.
#' @return numeric vector of roots.
#' @keywords internal
tsce_H_inverse <- function(z, lo, hi, k, tol = 1e-12, max_newton = 50L,
                           max_bisect = 200L) {
  t <- hi
  active <- seq_along(z)
  for (it in seq_len(max_newton)) {
    if (length(active) == 0L) break
    ka <- subset_k(k, active)
    ta <- t[active]
    f <- tsce_H(ta, ka) - z[active]
    h <- tsce_h(ta, ka)
    step <- f / h
    t_new <- pmax(ta - step, lo[active])
    bad <- !is.finite(t_new)
    t_new[bad] <- ta[bad]
    t[active] <- t_new
    done <- !bad & (abs(ta - t_new) <= tol * pmax(1, ta) | f <= 0)
    active <- active[!done]
  }
  if (length(active) > 0L) {
    # bisection on [lo, t], where H(t) >= z still holds
    a <- lo[active]
    b <- t[active]
    ka <- subset_k(k, active)
    za <- z[active]
    for (it in seq_len(max_bisect)) {
      m <- (a + b) / 2
      up <- tsce_H(m, ka) >= za
      b[up] <- m[up]
      a[!up] <- m[!up]
      if (all(b - a <= tol * pmax(1, b))) break
    }
    t[active] <- (a + b) / 2
  }
  t
}

#' Inverse of the cumulative onset intensity
#'
#' For each target \eqn{z} (cumulative intensity), finds \eqn{t} with
#' \eqn{\Lambda_{mal}(t) = z}: selects the smoking segment by comparing \eqn{z}
#' with \eqn{\Lambda_{mal}(t_s)} and \eqn{\Lambda_{mal}(t_q)}, forms the local
#' target \eqn{z' = z - \Lambda_{mal}(b_j) + \tilde H_{\psi_j}(b_j)}, and solves
#' \eqn{\tilde H_{\psi_j}(t) = z'} by [tsce_H_inverse()] on the segment, capped
#' at `args$shared$t_max`. Targets above \eqn{\Lambda_{mal}(t_{max})} (no event
#' by `t_max`) return `NA`; [nhppp::vdraw_cumulative_intensity()] does not pass
#' such targets.
#'
#' @param z numeric vector or matrix of targets, `NA` allowed; row `i` belongs
#'   to person `i`.
#' @param args list with `shared$t_max` and `row_args` (see [Lambda_mal()]).
#' @return numeric vector or matrix like `z`.
#' @keywords internal
Lambda_mal_inv <- function(z, args) {
  if (is.matrix(z)) {
    out <- z
    for (j in seq_len(ncol(z))) out[, j] <- Lambda_mal_inv(z[, j], args)
    return(out)
  }
  ra <- args$row_args
  t_max <- args$shared$t_max
  out <- rep(NA_real_, length(z))
  ok <- which(!is.na(z))
  if (length(ok) == 0L) {
    return(out)
  }
  sub <- list(row_args = ra[ok])
  zz <- z[ok]
  b1 <- pmin(ra$ts[ok], t_max)
  b2 <- pmin(ra$tq[ok], t_max)
  L1 <- Lambda_mal(b1, sub)
  L2 <- Lambda_mal(b2, sub)
  L_max <- Lambda_mal(rep(t_max, length(ok)), sub)
  k_ns <- list(B = ra$B_ns[ok], c = ra$c_ns[ok], K = ra$K_ns[ok], nmx = ra$nmx_ns[ok])
  k_s <- list(B = ra$B_s[ok], c = ra$c_s[ok], K = ra$K_s[ok], nmx = ra$nmx_s[ok])

  res <- rep(NA_real_, length(zz))
  seg1 <- zz <= L1
  seg2 <- !seg1 & zz <= L2
  beyond <- zz > L_max # no event in [0, t_max]
  seg3 <- !seg1 & !seg2 & !beyond
  if (any(seg1)) {
    i <- which(seg1)
    res[i] <- tsce_H_inverse(zz[i], rep(0, length(i)), b1[i], subset_k(k_ns, i))
  }
  if (any(seg2)) {
    i <- which(seg2)
    ks <- subset_k(k_s, i)
    zp <- zz[i] - L1[i] + tsce_H(b1[i], ks)
    res[i] <- tsce_H_inverse(zp, b1[i], b2[i], ks)
  }
  if (any(seg3)) {
    i <- which(seg3)
    kn <- subset_k(k_ns, i)
    zp <- zz[i] - L2[i] + tsce_H(b2[i], kn)
    res[i] <- tsce_H_inverse(zp, b2[i], rep(t_max, length(i)), kn)
  }
  out[ok] <- res
  out
}

#' Per-person arguments for [Lambda_mal()] and [Lambda_mal_inv()]
#'
#' @inheritParams generate_onset
#' @return data.table with `ts`, `tq` (`Inf` where absent) and the TSCE
#'   constants of the non-smoking (`_ns`) and smoking (`_s`) parameter sets.
#' @keywords internal
onset_row_args <- function(params, age_smoking_start, age_smoking_quit,
                           intensity) {
  k <- tsce_person_constants(params, intensity)
  ts <- age_smoking_start
  tq <- age_smoking_quit
  ts[is.na(ts)] <- Inf
  tq[is.na(tq)] <- Inf
  stopifnot(all(ts <= tq), all(ts >= 0))
  data.table(
    ts = ts, tq = tq,
    B_ns = k$ns$B, c_ns = k$ns$c, K_ns = k$ns$K, nmx_ns = k$ns$nmx,
    B_s = k$s$B, c_s = k$s$c, K_s = k$s$K, nmx_s = k$s$nmx
  )
}

#' Generate age at onset of the first malignant cell
#'
#' Samples \eqn{T_{mal}}, the first arrival of the TSCE onset process from
#' birth, with [nhppp::vdraw_cumulative_intensity()] on \eqn{[0, M]}, using the
#' closed-form cumulative intensity [Lambda_mal()] and its vectorized Newton
#' inverse [Lambda_mal_inv()]. This is the inverse-CDF map of the archival
#' `t_mal`. Onset is drawn from birth and is not conditioned on onset after
#' study entry (archival).
#'
#' @param params data.frame with one row per person holding the TSCE columns of
#'   [deterministic_params]; see [tsce_person_constants()].
#' @param age_smoking_start age at starting smoking (years); `NA` for never
#'   smokers.
#' @param age_smoking_quit age at quitting (years); `NA` for never and current
#'   smokers.
#' @param intensity cigarettes per day while smoking.
#' @param max_age \eqn{M}, the upper end of the sampling interval; default
#'   `getOption("milc.max_age")`. The archival code searches up to 150.
#' @return numeric vector, age at onset; `NA` when there is no onset by
#'   `max_age`.
#' @export
generate_onset <- function(params, age_smoking_start, age_smoking_quit,
                           intensity, max_age = getOption("milc.max_age")) {
  n <- nrow(params)
  stopifnot(
    length(age_smoking_start) == n, length(age_smoking_quit) == n,
    length(intensity) == n, is.numeric(max_age), length(max_age) == 1L
  )
  if (n == 0L) {
    return(numeric(0))
  }
  ra <- onset_row_args(params, age_smoking_start, age_smoking_quit, intensity)
  nhppp::vdraw_cumulative_intensity(
    Lambda = Lambda_mal,
    Lambda_inv = Lambda_mal_inv,
    t_min = 0,
    t_max = rep(max_age, n),
    Lambda_args = list(shared = list(t_max = max_age), row_args = ra),
    report_first_K = 1L
  )[, 1L]
}
