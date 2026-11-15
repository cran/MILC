# Tests for the TSCE onset kernel, its inverse, and generate_onset().

archival_onset <- function() {
  dir <- testthat::test_path("..", "..", "archival_code", "R")
  if (!dir.exists(dir)) {
    return(NULL)
  }
  e <- new.env()
  for (f in c("HT_mal.R", "ht_mal_int.R", "t_mal.R")) {
    sys.source(file.path(dir, f), envir = e)
  }
  e
}

# Profiles: never (male, female), current, former, heavy former smoker.
profiles <- data.frame(
  sex = c("male", "female", "male", "female", "male", "female"),
  ts = c(NA, NA, 16, 18, 15, 20),
  tq = c(NA, NA, NA, 50, 40, NA),
  q = c(0, 0, 20, 30, 201, 5)
)

profile_args <- function(rows = seq_len(nrow(profiles)), t_max = 100) {
  pr <- profiles[rows, ]
  pp <- deterministic_params[match(pr$sex, deterministic_params$sex)]
  list(
    shared = list(t_max = t_max),
    row_args = onset_row_args(pp, pr$ts, pr$tq, pr$q)
  )
}

archival_Lambda <- function(e, i, t) {
  g <- profiles$sex[i]
  d <- profiles$q[i]
  ts <- profiles$ts[i]
  tq <- profiles$tq[i]
  if (is.na(ts) || t <= ts) {
    return(e$ht_mal_int(0, t, g, d, "N"))
  }
  if (is.na(tq) || t <= tq) {
    return(e$ht_mal_int(0, ts, g, d, "N") + e$ht_mal_int(ts, t, g, d, "Y"))
  }
  e$ht_mal_int(0, ts, g, d, "N") + e$ht_mal_int(ts, tq, g, d, "Y") +
    e$ht_mal_int(tq, t, g, d, "N")
}

test_that("stable helpers agree with the direct formulas", {
  u <- c(-0.5, -0.05, -1e-4, 1e-8, 1e-3, 0.09, 0.2, 3)
  big <- abs(u) > 0.01
  expect_equal(expm1_minus_x(u)[big], (expm1(u) - u)[big], tolerance = 1e-12)
  expect_equal(expm1_minus_x(1e-8), 5e-17, tolerance = 1e-10)
  x <- c(1e-9, 1e-4, 0.005, 0.02, 1, 1e6)
  expect_equal(log1p_minus_x(x)[x > 0.01], (log1p(x) - x)[x > 0.01], tolerance = 1e-12)
  expect_equal(log1p_minus_x(1e-9), -5e-19, tolerance = 1e-8)
  expect_true(all(is.na(expm1_minus_x(c(NA, 1))[1])))
})

test_that("tsce_h is the derivative of tsce_H", {
  pp <- deterministic_params
  k <- tsce_person_constants(pp, c(20, 30))$s
  for (t in c(1, 20, 60, 95)) {
    eps <- 1e-4 * t
    num <- (tsce_H(rep(t + eps, 2), k) - tsce_H(rep(t - eps, 2), k)) / (2 * eps)
    expect_equal(tsce_h(rep(t, 2), k), num, tolerance = 1e-6)
  }
  expect_identical(tsce_H(c(0, 0), k), c(0, 0))
  expect_identical(tsce_h(c(0, 0), k), c(0, 0))
})

test_that("Lambda_mal matches the archival ht_mal_int over the smoking history", {
  e <- archival_onset()
  skip_if(is.null(e), "archival code not available")
  a <- profile_args()
  for (t in c(1, 10, 17, 30, 45, 60, 90)) {
    ours <- Lambda_mal(rep(t, nrow(profiles)), a)
    arch <- vapply(seq_len(nrow(profiles)), function(i) archival_Lambda(e, i, t), 1)
    expect_equal(ours, arch, tolerance = 1e-9)
  }
  # matrix input is evaluated column by column
  m <- matrix(c(rep(0, 6), rep(60, 6)), ncol = 2)
  expect_equal(Lambda_mal(m, a)[, 2], Lambda_mal(rep(60, 6), a))
  expect_identical(Lambda_mal(m, a)[, 1], rep(0, 6))
})

test_that("Lambda_mal_inv inverts Lambda_mal in every segment", {
  a <- profile_args()
  n <- nrow(profiles)
  ra <- a$row_args
  L_max <- Lambda_mal(rep(100, n), a)
  for (frac in c(1e-12, 1e-6, 0.01, 0.3, 0.7, 0.999)) {
    z <- frac * L_max
    t <- Lambda_mal_inv(z, a)
    expect_true(all(t >= 0 & t <= 100))
    expect_equal(Lambda_mal(t, a), z, tolerance = 1e-10)
  }
  # targets exactly at the segment boundaries
  z_ts <- Lambda_mal(pmin(ra$ts, 100), a)
  z_tq <- Lambda_mal(pmin(ra$tq, 100), a)
  expect_equal(Lambda_mal_inv(z_ts, a), pmin(ra$ts, 100), tolerance = 1e-10)
  expect_equal(Lambda_mal_inv(z_tq, a), pmin(ra$tq, 100), tolerance = 1e-10)
  # beyond t_max: no event; NA passes through
  expect_true(all(is.na(Lambda_mal_inv(L_max * 1.01, a))))
  expect_true(all(is.na(Lambda_mal_inv(rep(NA_real_, n), a))))
})

test_that("matched uniforms reproduce the archival t_mal", {
  e <- archival_onset()
  skip_if(is.null(e), "archival code not available")
  set.seed(20261121)
  U <- runif(150)
  for (i in seq_len(nrow(profiles))) {
    arch <- vapply(U, function(u) {
      e$t_mal(u, profiles$sex[i], profiles$ts[i], profiles$tq[i], profiles$q[i])
    }, 1)
    a <- profile_args(rep(i, length(U)), t_max = 150)
    ours <- Lambda_mal_inv(-log(U), a)
    # archival returns 150 when the root is beyond 150; we return NA
    expect_identical(is.na(ours), arch >= 150)
    ok <- !is.na(ours)
    # uniroot's default tolerance is about 1.2e-4
    expect_lt(max(abs(ours[ok] - arch[ok])), 1e-3)
  }
})

test_that("generate_onset draws from 1 - exp(-Lambda_mal) (Monte Carlo)", {
  set.seed(20261122)
  n <- 1e5L
  for (i in c(3L, 4L)) {
    a <- profile_args(rep(i, n))
    x <- generate_onset(
      deterministic_params[match(rep(profiles$sex[i], n), sex)],
      rep(profiles$ts[i], n), rep(profiles$tq[i], n), rep(profiles$q[i], n)
    )
    expect_true(all(is.na(x) | (x > 0 & x <= 100)))
    t <- c(30, 50, 70, 90, 100)
    target <- 1 - exp(-Lambda_mal(t, profile_args(rep(i, length(t)))))
    emp <- vapply(t, function(s) mean(!is.na(x) & x <= s), 1)
    se <- sqrt(target * (1 - target) / n)
    expect_true(all(abs(emp - target) <= 5 * se + 1e-12))
  }
})

test_that("generate_onset input handling", {
  P <- deterministic_params
  expect_identical(
    generate_onset(P[0], numeric(0), numeric(0), numeric(0)), numeric(0)
  )
  expect_error(generate_onset(P[1], 30, 20, 10)) # quit before start
  expect_error(generate_onset(P[1:2], 30, NA, 10)) # length mismatch
  set.seed(1)
  x <- generate_onset(P[c(1, 2)], c(NA, 18), c(NA, NA), c(0, 20))
  expect_length(x, 2)
})
