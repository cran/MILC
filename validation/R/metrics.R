# Effect-size comparisons for the validation harness.
#
# Three kinds of statistics:
#   - one-sample: sup over the table knots of |F_hat(t_k) - F_ref(t_k)|, with
#     F_ref the exact archival law of a mortality stratum;
# and two-sample, fastMILC minus archival:
#   - sup-norm distance between two (possibly defective) distribution
#     functions, D = sup_t |F_fast(t) - F_arch(t)|, where
#     F(t) = (1/n) #{i: X_i <= t}, X_i = NA for "no event" (defective CDF,
#     i.e. a cumulative incidence function without censoring);
#   - difference of two proportions.
# Each comes with an approximate 95% Monte Carlo bound b on |estimate - truth|:
#   - one-sample: b = 1.358 / sqrt(n) (Kolmogorov; conservative over knots);
#   - sup-norm: b = 1.358 sqrt(1/n_fast + 1/n_arch) (asymptotic Kolmogorov
#     quantile; exact under equality of the laws, approximate near it);
#   - proportion: b = 1.96 se (Agresti-Caffo adjusted).
# Verdict against a tolerance tol on the true effect size:
#   "equivalent"   if |est| + b <= tol
#   "different"    if |est| - b >  tol
#   "inconclusive" otherwise (Monte Carlo error too large to decide).
# Bounds are per row, not adjusted for multiplicity. Rows with role
# "descriptive" (comparisons with the unmodified archival samplers, which carry
# the D8 lookup lag) report the estimate and bound without a verdict.

TOL_CDF <- 0.005 # sup-norm
TOL_PROP <- 0.005 # proportions

verdict <- function(est, bound, tol) {
  a <- abs(est)
  data.table::fifelse(a + bound <= tol, "equivalent",
    data.table::fifelse(a - bound > tol, "different", "inconclusive")
  )
}

# x1, x2: values (NA = no event); denominators are the full lengths.
# Times are rounded to 1e-8 so that point masses computed in floating point
# on the two sides (e.g. T_diagn + h - T_diagn) coincide.
sup_diff <- function(x1, x2) {
  n1 <- length(x1)
  n2 <- length(x2)
  x1 <- round(x1, 8)
  x2 <- round(x2, 8)
  s1 <- sort(x1[!is.na(x1)])
  s2 <- sort(x2[!is.na(x2)])
  grid <- sort(unique(c(s1, s2)))
  if (!length(grid)) return(0)
  max(abs(findInterval(grid, s1) / n1 - findInterval(grid, s2) / n2))
}

cdf_row <- function(check, group, x_fast, x_arch, tol = TOL_CDF, role = "test") {
  n1 <- length(x_fast)
  n2 <- length(x_arch)
  est <- sup_diff(x_fast, x_arch)
  b <- 1.358 * sqrt(1 / n1 + 1 / n2)
  data.table::data.table(
    role = role, check = check, group = group, statistic = "sup|dF|",
    n_fast = n1, n_arch = n2, value_fast = NA_real_, value_arch = NA_real_,
    estimate = est, bound = b, tol = tol,
    verdict = if (role == "test") verdict(est, b, tol) else "descriptive"
  )
}

prop_row <- function(check, group, y_fast, y_arch, tol = TOL_PROP, role = "test") {
  n1 <- length(y_fast)
  n2 <- length(y_arch)
  p1 <- mean(y_fast)
  p2 <- mean(y_arch)
  est <- p1 - p2
  # Agresti-Caffo: add one success and one failure per side for the se
  q1 <- (sum(y_fast) + 1) / (n1 + 2)
  q2 <- (sum(y_arch) + 1) / (n2 + 2)
  b <- 1.96 * sqrt(q1 * (1 - q1) / (n1 + 2) + q2 * (1 - q2) / (n2 + 2))
  data.table::data.table(
    role = role, check = check, group = group, statistic = "dp",
    n_fast = n1, n_arch = n2, value_fast = p1, value_arch = p2,
    estimate = est, bound = b, tol = tol,
    verdict = if (role == "test") verdict(est, b, tol) else "descriptive"
  )
}

# One-sample: draws x (time since entry or diagnosis; NA = none) against the
# exact reference CDF F_ref at knots t (t > 0). `n_na`: draws that returned NA
# (archival lookup), reported separately.
knot_row <- function(check, group, x, t, F_ref, sampler, role = "test",
                     tol = TOL_CDF) {
  n <- length(x)
  s <- sort(x[!is.na(x)])
  est <- max(abs(findInterval(t, s) / n - F_ref))
  b <- 1.358 / sqrt(n)
  data.table::data.table(
    role = role, check = check, group = group, sampler = sampler,
    statistic = "max_k|F_hat-F_ref|", n = n, na_rate = mean(is.na(x)),
    estimate = est, bound = b, tol = tol,
    verdict = if (role == "test") verdict(est, b, tol) else "descriptive"
  )
}

# Assembly of a latent history (as in Population$.assembleHistory): death is
# the first latent death time (lung cancer if T_dl <= T_do); a death after M
# means alive at M; progression events after the end of follow-up are not
# observed.
assemble <- function(d, M) {
  d[, age_dead := pmin(T_do, T_dl, na.rm = TRUE)]
  d[, cause := data.table::fifelse(!is.na(T_dl) & T_dl <= T_do, "lung_cancer", "other_causes")]
  d[age_dead > M, `:=`(age_dead = NA_real_, cause = NA_character_)]
  end <- d[, data.table::fifelse(is.na(age_dead), M, age_dead)]
  for (v in c("T_mal", "T_reg", "T_dist", "T_diagn")) {
    data.table::set(d, which(!(d[[v]] <= end)), v, NA_real_)
  }
  d[is.na(T_diagn), `:=`(D_diagn = NA_real_, stage = NA_character_)]
  d[]
}

intensity_group <- function(q) {
  as.character(cut(q, c(-Inf, 0, 10, 20, Inf), c("0", "1-10", "11-20", "21+")))
}

# Subgroups for end-to-end checks: all; sex x status; current smokers by
# sex x intensity group. Returns a list of (label, logical index) per side.
subgroups <- function(d) {
  g <- list(all = rep(TRUE, nrow(d)))
  for (s in c("female", "male")) for (st in c("former", "current")) {
    g[[paste(s, st)]] <- d$sex == s & d$status == st
  }
  for (s in c("female", "male")) for (q in c("1-10", "11-20", "21+")) {
    g[[paste(s, "current", q)]] <- d$sex == s & d$status == "current" & d$intgrp == q
  }
  g
}

# End-to-end endpoints. `f`, `a`: assembled histories with columns sex,
# status, intgrp, entry, T_mal, T_diagn, D_diagn, stage, age_dead, cause,
# T_dl (latent).
# Subgroups with fewer than `min_n` people on either side are skipped.
end_to_end_metrics <- function(f, a, role = "test", horizons = c(5, 10, 20),
                               min_n = 1000L) {
  gf <- subgroups(f)
  ga <- subgroups(a)
  rows <- list()
  add <- function(x) rows[[length(rows) + 1L]] <<- x
  cdf_r <- function(...) cdf_row(..., role = role)
  prop_r <- function(...) prop_row(..., role = role)
  for (g in names(gf)) {
    F <- f[gf[[g]]]
    A <- a[ga[[g]]]
    if (min(nrow(F), nrow(A)) < min_n) next
    tdx <- function(d) d$T_diagn - d$entry
    tdl <- function(d) data.table::fifelse(d$cause == "lung_cancer", d$age_dead - d$entry, NA_real_)
    tdo <- function(d) data.table::fifelse(d$cause == "other_causes", d$age_dead - d$entry, NA_real_)
    add(cdf_r("onset (observed), by age", g, F$T_mal, A$T_mal, TOL_CDF))
    add(cdf_r("diagnosis, by time since entry", g, tdx(F), tdx(A), TOL_CDF))
    add(cdf_r("lung-cancer death, by time since entry", g, tdl(F), tdl(A), TOL_CDF))
    add(cdf_r("other-cause death, by time since entry", g, tdo(F), tdo(A), TOL_CDF))
    for (h in horizons) {
      add(prop_r(sprintf("diagnosed within %d y", h), g,
                   !is.na(tdx(F)) & tdx(F) <= h, !is.na(tdx(A)) & tdx(A) <= h))
      add(prop_r(sprintf("lung-cancer death within %d y", h), g,
                   !is.na(tdl(F)) & tdl(F) <= h, !is.na(tdl(A)) & tdl(A) <= h))
      add(prop_r(sprintf("other-cause death within %d y", h), g,
                   !is.na(tdo(F)) & tdo(F) <= h, !is.na(tdo(A)) & tdo(A) <= h))
    }
    add(prop_r("diagnosed, lifetime", g, !is.na(F$T_diagn), !is.na(A$T_diagn)))
    add(prop_r("lung-cancer death, lifetime", g,
                 F$cause %in% "lung_cancer", A$cause %in% "lung_cancer"))
    Fd <- F[!is.na(T_diagn)]
    Ad <- A[!is.na(T_diagn)]
    for (s in c("localized", "regional", "distant")) {
      add(prop_r(paste("stage share:", s), g, Fd$stage == s, Ad$stage == s))
    }
    add(cdf_r("diameter at diagnosis", g, Fd$D_diagn, Ad$D_diagn, TOL_CDF))
    for (s in c("localized", "regional", "distant")) {
      add(cdf_r(paste("latent survival after diagnosis:", s), g,
                  Fd[stage == s, T_dl - T_diagn], Ad[stage == s, T_dl - T_diagn],
                  TOL_CDF))
    }
  }
  data.table::rbindlist(rows)
}

# Cumulative incidence curves (no censoring) on a time grid, for figures.
cif_curves <- function(d, side, grid = seq(0, 50, by = 0.25)) {
  out <- list()
  for (g in c("all", "former", "current")) {
    D <- if (g == "all") d else d[status == g]
    n <- nrow(D)
    for (ev in c("diagnosis", "lung_cancer", "other_causes")) {
      t <- if (ev == "diagnosis") D$T_diagn - D$entry else
        data.table::fifelse(D$cause %in% ev, D$age_dead - D$entry, NA_real_)
      s <- sort(t[!is.na(t)])
      out[[length(out) + 1L]] <- data.table::data.table(
        side = side, group = g, event = ev, time = grid,
        cif = findInterval(grid, s) / n
      )
    }
  }
  data.table::rbindlist(out)
}
