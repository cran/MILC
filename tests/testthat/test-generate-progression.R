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

test_that("progression ages are ordered and consistent with stage", {
  set.seed(1)
  n <- 2e4L
  onset <- c(rep(NA, 10), runif(n - 10, 20, 80))
  x <- generate_progression(onset, prog_params(n))
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

test_that("distribution matches the archival t_prog (swap included)", {
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
  ours <- generate_progression(rep(0, n), prog_params(n))
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
