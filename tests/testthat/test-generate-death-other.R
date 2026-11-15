# Tests for generate_death_other() and its helpers.

archival_r <- function(file) {
  testthat::test_path("..", "..", "archival_code", "R", file)
}

# CDF of the latent duration implied by a rate row, at durations `t` within the
# tabulated window (exact; test-mortality-data.R checks it against the tables).
rate_cdf <- function(rates, row, t) {
  lam <- rates$lambda_matrix[row, ]
  tb <- rates$time_breaks
  Lambda_at <- function(x) {
    j <- findInterval(x, tb, rightmost.closed = TRUE)
    cum <- c(0, cumsum(lam * diff(tb)))
    cum[j] + lam[j] * (x - tb[j])
  }
  1 - exp(-vapply(t, Lambda_at, numeric(1)))
}

expect_cdf <- function(x, t, target, n_se = 5) {
  n <- length(x)
  emp <- vapply(t, function(s) mean(x <= s), numeric(1))
  se <- sqrt(target * (1 - target) / n)
  expect_true(
    all(abs(emp - target) <= n_se * se + 1e-12),
    info = paste("emp:", paste(signif(emp, 4), collapse = " "),
                 "target:", paste(signif(target, 4), collapse = " "))
  )
}

draw_one_stratum <- function(n, age, sex, status, q, ...) {
  generate_death_other(
    rep(age, n), rep(sex, n), rep(status, n), rep(q, n), ...
  )
}

test_that("age and intensity groups match the archival scalar functions", {
  skip_if_not(file.exists(archival_r("age_grp.R")), "archival code not available")
  env <- new.env()
  sys.source(archival_r("age_grp.R"), envir = env)
  sys.source(archival_r("d_grp.R"), envir = env)
  ages <- c(0, 39.999, 40, 44.999, 45, 62.5, 84.999, 85, 100, 119.9, 120, 121)
  expect_identical(
    milc_age_group(ages),
    suppressWarnings(vapply(ages, env$age_grp, ""))
  )
  q <- c(0, 0.5, 1, 10.99, 11, 20.99, 21, 60)
  expect_identical(milc_intensity_group(q), vapply(q, env$d_grp, ""))
})

test_that("stratum lookup maps labels and rejects unknown strata", {
  s <- rates_death_other$strata
  expect_identical(
    death_other_stratum(41, "Male", "non_smoker", 0),
    s[smoking_status == "never" & agegrp == "40-44" & sex == "male", stratum_row]
  )
  expect_identical(
    death_other_stratum(41, "female", "smoker", 15),
    s[smoking_status == "current" & agegrp == "40-44" & sex == "female" &
        intensity == "11-20", stratum_row]
  )
  # intensity is ignored for never and former smokers
  expect_identical(
    death_other_stratum(c(50, 50), c("male", "male"), c("former", "former"), c(0, 30)),
    rep(s[smoking_status == "former" & agegrp == "50-54" & sex == "male",
          stratum_row], 2L)
  )
  expect_error(death_other_stratum(50, "male", "occasional", 0), "Unknown")
  # a current smoker with fewer than 1 cigarette/day has no table (archival)
  expect_error(death_other_stratum(50, "male", "current", 0), "No other-cause")
  # age 120 maps to NA (archival), which has no other-cause table
  expect_error(death_other_stratum(120, "male", "never", 0), "No other-cause")
})

test_that("corrected tail: window reproduces the table, tail is Unif(a+h, M)", {
  M <- getOption("milc.max_age")
  set.seed(20261115)
  n <- 2e5L
  age <- 70
  row <- death_other_stratum(age, "male", "never", 0)
  h <- rates_death_other$strata$horizon[row]
  x <- draw_one_stratum(n, age, "male", "never", 0, tail = "corrected") - age
  t <- c(0.5, 1, 2, 4, 8, 12, 16)
  expect_cdf(x, t, rate_cdf(rates_death_other, row, t))
  tail <- x[x > h]
  expect_true(all(tail <= M - age))
  F_max <- rates_death_other$strata$F_max[row]
  expect_cdf(
    x, h + c(2, 6, 10),
    F_max + (1 - F_max) * c(2, 6, 10) / (M - age - h)
  )
})

test_that("archival tail: mixture of the table and Unif(a, M - h)", {
  M <- getOption("milc.max_age")
  set.seed(20261116)
  n <- 2e5L
  age <- 42
  row <- death_other_stratum(age, "female", "never", 0)
  h <- rates_death_other$strata$horizon[row]
  F_max <- rates_death_other$strata$F_max[row]
  x <- draw_one_stratum(n, age, "female", "never", 0) - age
  expect_true(all(x > 0 & x <= M - h - age))
  t <- c(1, 4, 8, 16, 30, 50)
  target <- (1 - F_max) * pmin(1, t / (M - age - h))
  target[t <= h] <- target[t <= h] + rate_cdf(rates_death_other, row, t[t <= h])
  target[t > h] <- target[t > h] + F_max
  expect_cdf(x, t, target)
})

test_that("archival tail reproduces death before entry for late entry", {
  set.seed(1)
  M <- getOption("milc.max_age")
  x <- draw_one_stratum(2e4L, 95, "male", "never", 0)
  # M - 16 < 95: archival tail draws fall in (M - 16, 95)
  expect_true(any(x < 95))
  expect_true(all(x > M - 16 - 1e-12))
  expect_warning(
    draw_one_stratum(100L, 95, "male", "never", 0, tail = "corrected"),
    "reversed"
  )
})

test_that("terminal-certain stratum: no tail, uniform in the terminal cell", {
  set.seed(20261117)
  n <- 2e5L
  row <- death_other_stratum(87, "female", "current", 25)
  st <- rates_death_other$strata[row]
  expect_true(st$terminal_certain)
  x <- draw_one_stratum(n, 87, "female", "current", 25) - 87
  expect_true(all(x >= 0 & x <= st$horizon))
  t <- c(0.5, 1, st$terminal_lo)
  expect_cdf(x, t, rate_cdf(rates_death_other, row, t))
  last <- x[x > st$terminal_lo]
  w <- st$horizon - st$terminal_lo
  expect_cdf(last, st$terminal_lo + w * c(0.25, 0.5, 0.75), c(0.25, 0.5, 0.75))
})

test_that("mixed population: chunking and stratum order do not change the law", {
  pop <- data.table::data.table(
    age = c(rep(45, 3000), rep(70, 3000), rep(87, 3000)),
    sex = rep(c("male", "female", "female"), each = 3000),
    status = rep(c("former", "never", "current"), each = 3000),
    q = rep(c(0, 0, 25), each = 3000)
  )
  set.seed(5)
  a <- pop[, generate_death_other(age, sex, status, q, chunk_size = 1e5)]
  set.seed(5)
  b <- pop[, generate_death_other(age, sex, status, q, chunk_size = 1e5)]
  expect_identical(a, b)
  set.seed(5)
  small <- pop[, generate_death_other(age, sex, status, q, chunk_size = 7)]
  expect_length(small, nrow(pop))
  expect_false(anyNA(small))
  # same group means within Monte Carlo error
  ma <- tapply(a, pop$age, mean)
  ms <- tapply(small, pop$age, mean)
  sds <- tapply(a, pop$age, sd) * sqrt(2 / 3000)
  expect_true(all(abs(ma - ms) <= 5 * sds))
})

test_that("input checks", {
  expect_error(generate_death_other(50, "male", "never", 0, tail = "other"))
  expect_error(generate_death_other(NA_real_, "male", "never", 0))
  expect_error(generate_death_other(c(50, 60), "male", "never", 0))
  expect_identical(
    generate_death_other(numeric(0), character(0), character(0), numeric(0)),
    numeric(0)
  )
})
