# Tests for generate_progression() and its helpers.

archival_prog <- function() {
  dir <- testthat::test_path("..", "..", "archival_code", "R")
  if (!dir.exists(dir)) {
    return(NULL)
  }
  e <- new.env()
  for (f in c("t_vol.R", "t_prog.R")) sys.source(file.path(dir, f), envir = e)
  e
}

prog_params <- function(n) {
  deterministic_params[rep(1L, n)]
}

test_that("milc_time_to_volume matches the archival t_vol", {
  e <- archival_prog()
  skip_if(is.null(e), "archival code not available")
  v <- c(1e-6, 1, 1e3, 5e5, 1e6)
  m <- c(3.4e-4, 4.2e-4, 1e-3, 3.4e-4, 2e-4)
  expect_equal(
    milc_time_to_volume(v, m),
    mapply(e$t_vol, v, m),
    tolerance = 1e-12
  )
})

test_that("truncated volumes stay in (V0, Vmax) and follow the truncated law", {
  set.seed(20261123)
  n <- 1e5L
  v <- draw_truncated_volume(rep(5.62, n), rep(5.62, n))
  expect_true(all(v > MILC_V0 & v < MILC_VMAX))
  lo <- pnorm((log(MILC_V0 / 1000) - 5.62) / 5.62)
  hi <- pnorm((log(MILC_VMAX / 1000) - 5.62) / 5.62)
  q <- c(1e-3, 1, 1e3, 1e5)
  target <- (pnorm((log(q / 1000) - 5.62) / 5.62) - lo) / (hi - lo)
  emp <- vapply(q, function(s) mean(v <= s), 1)
  expect_true(all(abs(emp - target) <= 5 * sqrt(target * (1 - target) / n) + 1e-12))
})

for (method in c("rejection", "conditional")) test_that(paste("progression ages are ordered and consistent with stage:", method), {
  set.seed(1)
  n <- 2e4L
  onset <- c(rep(NA, 10), runif(n - 10, 20, 80))
  x <- generate_progression(onset, prog_params(n), method = method)
  expect_true(all(is.na(x[1:10])))
  y <- x[-(1:10)]
  o <- onset[-(1:10)]
  expect_false(anyNA(y))
  expect_true(all(y$age_at_regional_disease > o))
  expect_true(all(y$age_at_distal_disease > y$age_at_regional_disease))
  expect_true(all(y$age_at_diagnosis > o))
  expect_true(all(y$diameter_at_diagnosis > 0.01 & y$diameter_at_diagnosis < 130))
  stage <- with(y, fifelse(
    age_at_diagnosis < age_at_regional_disease, "localized",
    fifelse(age_at_diagnosis < age_at_distal_disease, "regional", "distant")
  ))
  expect_identical(y$stage_at_diagnosis, stage)
  expect_setequal(unique(y$stage_at_diagnosis), c("localized", "regional", "distant"))
})

for (method in c("rejection", "conditional")) test_that(paste("distribution matches the archival t_prog (swap included):", method), {
  e <- archival_prog()
  skip_if(is.null(e), "archival code not available")
  p <- deterministic_params[1]
  cdiagn <- c(p$volume_diagnosis_meanlog, p$volume_diagnosis_sdlog)
  creg <- c(p$volume_regional_meanlog, p$volume_regional_sdlog)
  cdist <- c(p$volume_distant_meanlog, p$volume_distant_sdlog)
  n <- 4000L
  set.seed(20261124)
  arch <- data.table::rbindlist(lapply(seq_len(n), function(i) {
    r <- e$t_prog(1, p$tumor_growth_rate, cdiagn, creg, cdist)
    data.table::data.table(
      reg = r$Treg, dist = r$Tdist, diagn = r$Tdiagn, d = r$Ddiagn, stage = r$stage
    )
  }))
  set.seed(20261125)
  ours <- generate_progression(rep(0, n), prog_params(n), method = method)
  ks_p <- c(
    stats::ks.test(arch$reg, ours$age_at_regional_disease)$p.value,
    stats::ks.test(arch$dist, ours$age_at_distal_disease)$p.value,
    stats::ks.test(arch$diagn, ours$age_at_diagnosis)$p.value,
    stats::ks.test(arch$d, ours$diameter_at_diagnosis)$p.value
  )
  expect_true(all(ks_p > 1e-3), info = paste(signif(ks_p, 3), collapse = " "))
  tab <- rbind(
    table(factor(arch$stage, c("localized", "regional", "distant"))),
    table(factor(ours$stage_at_diagnosis, c("localized", "regional", "distant")))
  )
  expect_gt(stats::chisq.test(tab)$p.value, 1e-3)
})

# Conditional (common-random-numbers) sampler ------------------------------

test_that("Gauss-Legendre rule integrates polynomials exactly", {
  gl <- gauss_legendre(10L)
  expect_equal(sum(gl$w), 2, tolerance = 1e-14)
  expect_equal(sum(gl$w * gl$x^18), 2 / 19, tolerance = 1e-13)
  expect_equal(sum(gl$w * gl$x^7), 0, tolerance = 1e-14)
})

test_that("qtruncnorm is the truncated-normal quantile in both tails", {
  u <- c(1e-6, 0.1, 0.5, 0.9, 1 - 1e-6)
  for (w in list(c(-30, -25), c(-1, 2), c(9, 12))) {
    q <- qtruncnorm(u, 0, 1, w[1], w[2])
    expect_true(all(q >= w[1] & q <= w[2]))
    cdf <- pnorm_interval(w[1], q, 0, 1) / pnorm_interval(w[1], w[2], 0, 1)
    expect_equal(cdf, u, tolerance = 1e-8)
  }
  expect_length(qtruncnorm(0.5, 1, 2, -3, c(0, 1, 2)), 3L)
})

test_that("marginal of the larger ordered normal: table and inverse", {
  a <- MILC_LOGZ_MIN
  b <- MILC_LOGZ_MAX
  tab <- ordered_marginal_table(2.16, 2.16, 5.62, 5.62)
  g <- function(x) dnorm(x, 2.16, 2.16) * (pnorm(x, 5.62, 5.62) - pnorm(a, 5.62, 5.62))
  Ib <- integrate(g, a, b, rel.tol = 1e-13)$value
  expect_equal(tab$I[length(tab$I)], Ib, tolerance = 1e-10)
  u <- c(1e-8, 1e-3, 0.2, 0.5, 0.8, 0.999, 1 - 1e-9)
  x <- ordered_marginal_q(u, tab)
  expect_true(all(diff(x) > 0))
  expect_true(all(x > a & x < b))
  back <- vapply(x, function(s) integrate(g, a, s, rel.tol = 1e-13)$value, 1) / Ib
  expect_equal(back, u, tolerance = 1e-10)
})

test_that("conditional draws reproduce the truncated joint law (vs rejection)", {
  a <- MILC_LOGZ_MIN
  b <- MILC_LOGZ_MAX
  set.seed(20261126)
  N <- 2e6L
  xr <- rnorm(N, 5.62, 5.62) # V_reg uses cdist (archival swap)
  xd <- rnorm(N, 2.16, 2.16) # V_dist uses creg
  ok <- xr > a & xr < xd & xd < b
  xr <- xr[ok]
  xd <- xd[ok]
  n <- 50000L
  y <- generate_progression(rep(0, n), prog_params(n), method = "conditional")
  m <- deterministic_params$tumor_growth_rate[1]
  p_reg <- suppressWarnings(stats::ks.test(
    y$age_at_regional_disease, milc_time_to_volume(1000 * exp(xr), m)
  )$p.value)
  p_dist <- suppressWarnings(stats::ks.test(
    y$age_at_distal_disease, milc_time_to_volume(1000 * exp(xd), m)
  )$p.value)
  expect_gt(p_reg, 1e-3)
  expect_gt(p_dist, 1e-3)
  # dependence: both below their medians
  target <- mean(xr < median(xr) & xd < median(xd))
  yr <- y$age_at_regional_disease
  yd <- y$age_at_distal_disease
  emp <- mean(yr < median(yr) & yd < median(yd))
  expect_lt(abs(emp - target), 5 * sqrt(target * (1 - target) / n))
})

test_that("conditional method: three uniforms per person, monotone in parameters", {
  set.seed(7)
  invisible(runif(15))
  after <- runif(1)
  set.seed(7)
  invisible(generate_progression(rep(50, 5), prog_params(5), method = "conditional"))
  expect_identical(runif(1), after)
  # common random numbers: a larger growth rate shortens every time
  P1 <- prog_params(200)
  P2 <- data.table::copy(P1)[, tumor_growth_rate := tumor_growth_rate * 1.1]
  set.seed(8)
  x1 <- generate_progression(rep(50, 200), P1, method = "conditional")
  set.seed(8)
  x2 <- generate_progression(rep(50, 200), P2, method = "conditional")
  expect_true(all(x2$age_at_diagnosis < x1$age_at_diagnosis))
  expect_identical(x1$stage_at_diagnosis, x2$stage_at_diagnosis)
})

test_that("method is validated", {
  expect_error(generate_progression(50, prog_params(1), method = "other"))
})
