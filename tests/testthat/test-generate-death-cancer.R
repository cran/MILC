# Tests for generate_death_cancer() and its helpers.

lung_row <- function(st_, ag_, sx_) {
  rates_death_lung$strata[
    stage == st_ & agegrp == ag_ & sex == sx_,
    stratum_row
  ]
}

# CDF of the latent duration implied by a rate row, at durations `t` within the
# tabulated window.
lung_rate_cdf <- function(row, t) {
  lam <- rates_death_lung$lambda_matrix[row, ]
  tb <- rates_death_lung$time_breaks
  cum <- c(0, cumsum(lam * diff(tb)))
  j <- findInterval(t, tb, rightmost.closed = TRUE)
  1 - exp(-(cum[j] + lam[j] * (t - tb[j])))
}

expect_cdf_lung <- function(x, t, target, n_se = 5) {
  n <- length(x)
  emp <- vapply(t, function(s) mean(x <= s), numeric(1))
  se <- sqrt(target * (1 - target) / n)
  expect_true(
    all(abs(emp - target) <= n_se * se + 1e-12),
    info = paste("emp:", paste(signif(emp, 4), collapse = " "),
                 "target:", paste(signif(target, 4), collapse = " "))
  )
}

test_that("stage labels map to the archival labels", {
  expect_identical(
    milc_stage(c("local_disease", "regional_disease", "distal_disease",
                 "localized", "regional", "distant")),
    c("localized", "regional", "distant", "localized", "regional", "distant")
  )
  expect_error(milc_stage("in situ"), "Unknown stage")
})

test_that("stratum lookup uses age group at diagnosis and the NA fallback", {
  expect_identical(
    death_cancer_stratum(62, "distal_disease", "Male"),
    lung_row("distant", "60-64", "male")
  )
  expect_identical(
    death_cancer_stratum(c(39, 90), c("localized", "regional"), c("female", "male")),
    c(lung_row("localized", "<40", "female"), lung_row("regional", ">85", "male"))
  )
  # age exactly 120 -> archival age_grp NA -> literal "NA" strata
  expect_identical(
    death_cancer_stratum(120, "distant", "male"),
    lung_row("distant", "NA", "male")
  )
  expect_error(death_cancer_stratum(120, "regional", "male"), "No lung-cancer")
})

test_that("within the window the draws reproduce the table; residual at horizon", {
  set.seed(20261118)
  n <- 2e5L
  row <- lung_row("distant", "60-64", "male")
  st <- rates_death_lung$strata[row]
  x <- generate_death_cancer(rep(62, n), rep("distant", n), rep("male", n)) - 62
  expect_true(all(x > 0 & x <= st$horizon + 1e-12))
  t <- c(1, 3, 6, 12, 36, 120) / 12
  expect_cdf_lung(x, t, lung_rate_cdf(row, t))
  # residual mass 1 - F_max is a point mass at the horizon
  at_h <- mean(abs(x - st$horizon) < 1e-9)
  se <- sqrt(st$F_max * (1 - st$F_max) / n)
  expect_lt(abs(at_h - (1 - st$F_max)), 5 * se)
})

test_that("the merged atom: first-cell mass equals F(tau_1)", {
  skip_if_not(
    dir.exists(testthat::test_path("..", "..", "archival_code", "data")),
    "archival data not available"
  )
  env <- new.env()
  load(testthat::test_path("..", "..", "archival_code", "data", "ci.lung.rda"),
       envir = env)
  cif <- env$ci.lung[["regional 70-74 female 1"]]
  k <- sort(unique(cif$time))
  F1 <- max(cif$est[cif$time == k[2]])
  set.seed(20261119)
  n <- 2e5L
  x <- generate_death_cancer(rep(72, n), rep("regional", n), rep("female", n)) - 72
  expect_cdf_lung(x, k[2] / 12, F1)
})

test_that("degenerate strata: certain terminal cell and empty table", {
  set.seed(20261120)
  n <- 5e4L
  st <- rates_death_lung$strata[lung_row("distant", "NA", "male")]
  x <- generate_death_cancer(rep(120, n), rep("distant", n), rep("male", n)) - 120
  expect_true(all(x >= st$terminal_lo & x <= st$horizon))
  expect_cdf_lung(x, st$horizon * c(0.25, 0.5, 0.75), c(0.25, 0.5, 0.75))
  y <- generate_death_cancer(rep(120, 10), rep("localized", 10), rep("female", 10))
  expect_identical(y, rep(120, 10))
})

test_that("undiagnosed individuals get NA; stage is ignored for them", {
  set.seed(1)
  x <- generate_death_cancer(
    c(NA, 65, NA, 50), c(NA, "regional_disease", "bogus", "local_disease"),
    rep("female", 4)
  )
  expect_true(is.na(x[1]) && is.na(x[3]))
  expect_true(all(x[c(2, 4)] > c(65, 50)))
  expect_identical(
    generate_death_cancer(c(NA_real_, NA_real_), c(NA, NA), c("male", "male")),
    c(NA_real_, NA_real_)
  )
  expect_error(generate_death_cancer(60, NA_character_, "male"), "stage")
})

test_that("reproducible given the seed and independent of chunk size in law", {
  n <- 6000L
  age <- rep(c(55, 68, 77), each = n / 3)
  stage <- rep(c("localized", "regional", "distant"), each = n / 3)
  sex <- rep(c("male", "female"), n / 2)
  set.seed(9)
  a <- generate_death_cancer(age, stage, sex)
  set.seed(9)
  b <- generate_death_cancer(age, stage, sex)
  expect_identical(a, b)
  set.seed(9)
  s <- generate_death_cancer(age, stage, sex, chunk_size = 11)
  expect_false(anyNA(s))
  ma <- tapply(a - age, stage, mean)
  ms <- tapply(s - age, stage, mean)
  se <- tapply(a - age, stage, sd) * sqrt(2 / (n / 3))
  expect_true(all(abs(ma - ms) <= 5 * se))
})
