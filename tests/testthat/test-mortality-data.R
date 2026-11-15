# Tests for the step-rate tables rates_death_other and rates_death_lung.
# The exactness tests compare against the archival CIF tables, which live in the
# source tree only (archival_code/ is .Rbuildignore'd), so they skip elsewhere.

archival_dir <- function() {
  testthat::test_path("..", "..", "archival_code", "data")
}

load_archival <- function(names) {
  env <- new.env()
  for (nm in names) {
    load(file.path(archival_dir(), paste0(nm, ".rda")), envir = env)
  }
  env
}

# CIF implied by one row of a rate table, evaluated at the table's knots.
implied_cif <- function(rates, row, knots_years) {
  lam <- rates$lambda_matrix[row, ]
  Lambda <- c(0, cumsum(lam * diff(rates$time_breaks)))
  idx <- match(round(knots_years, 10), round(rates$time_breaks, 10))
  stopifnot(!anyNA(idx))
  1 - exp(-Lambda[idx])
}

# Max absolute difference between table and implied CIF over a family.
max_knot_error <- function(rates, tables, time_per_year) {
  errs <- vapply(seq_along(tables), function(i) {
    cif <- tables[[i]]
    k <- unique(cif$time)
    Fk <- vapply(k, function(t) max(cif$est[cif$time == t]), numeric(1))
    st <- rates$strata[i]
    keep <- k > 0 # the atom at 0 is merged into the first cell
    if (st$terminal_certain) {
      keep <- keep & k / time_per_year <= st$terminal_lo + 1e-12
    }
    if (!any(keep)) {
      return(0)
    }
    est <- implied_cif(rates, i, k[keep] / time_per_year)
    max(abs(est - Fk[keep]))
  }, numeric(1))
  max(errs)
}

test_that("rate tables have the structure vdraw_sc_step expects", {
  for (rates in list(rates_death_other, rates_death_lung)) {
    tb <- rates$time_breaks
    lm <- rates$lambda_matrix
    expect_true(is.numeric(tb) && all(diff(tb) > 0) && tb[1] == 0)
    expect_true(is.matrix(lm))
    expect_identical(ncol(lm), length(tb) - 1L)
    expect_identical(nrow(lm), nrow(rates$strata))
    expect_true(all(is.finite(lm)) && all(lm >= 0))
    expect_identical(rates$strata$stratum_row, seq_len(nrow(lm)))
  }
  expect_identical(dim(rates_death_other$lambda_matrix), c(110L, 64L))
  expect_identical(dim(rates_death_lung$lambda_matrix), c(68L, 245L))
})

test_that("rates are zero beyond each stratum's horizon", {
  for (rates in list(rates_death_other, rates_death_lung)) {
    left <- rates$time_breaks[-length(rates$time_breaks)]
    beyond <- outer(rates$strata$horizon, left, `<=`)
    expect_true(all(rates$lambda_matrix[beyond] == 0))
  }
})

test_that("other-cause rates reproduce the archival CIFs at every knot", {
  skip_if_not(dir.exists(archival_dir()), "archival data not available")
  env <- load_archival(c("never.other", "former.other", "current.other"))
  tables <- c(env$never.other, env$former.other, env$current.other)
  expect_identical(names(tables), rownames(rates_death_other$lambda_matrix))
  expect_lt(max_knot_error(rates_death_other, tables, 1), 1e-12)
})

test_that("lung rates reproduce the archival CIFs at every positive knot", {
  skip_if_not(dir.exists(archival_dir()), "archival data not available")
  env <- load_archival("ci.lung")
  expect_identical(names(env$ci.lung), rownames(rates_death_lung$lambda_matrix))
  expect_lt(max_knot_error(rates_death_lung, env$ci.lung, 12), 1e-12)
})

test_that("degenerate strata are flagged", {
  oth <- rates_death_other$strata
  lng <- rates_death_lung$strata
  expect_identical(
    oth[terminal_certain == TRUE, paste(smoking_status, agegrp, sex, intensity)],
    "current >85 female 21+"
  )
  expect_identical(
    lng[terminal_certain == TRUE, paste(stage, agegrp, sex)],
    "distant NA male"
  )
  expect_identical(lng[F_max == 0, paste(stage, agegrp, sex)], "localized NA female")
  expect_true(all(oth$F0 == 0))
  expect_true(all(lng$terminal_lo[lng$terminal_certain] < lng$horizon[lng$terminal_certain]))
})

test_that("vdraw_sc_step first arrivals match the tables (Monte Carlo)", {
  set.seed(20261115)
  n <- 200000L
  check_row <- function(rates, row) {
    x <- nhppp::vdraw_sc_step(
      lambda_matrix = rates$lambda_matrix[rep(row, n), , drop = FALSE],
      time_breaks = rates$time_breaks,
      report_first_K = 1L
    )[, 1]
    lam <- rates$lambda_matrix[row, ]
    Lambda <- c(0, cumsum(lam * diff(rates$time_breaks)))
    probe <- unique(round(seq(2, length(rates$time_breaks), length.out = 6)))
    target <- 1 - exp(-Lambda[probe])
    emp <- vapply(
      rates$time_breaks[probe],
      function(t) mean(!is.na(x) & x <= t),
      numeric(1)
    )
    se <- sqrt(target * (1 - target) / n)
    expect_true(all(abs(emp - target) <= 5 * se + 1e-12))
  }
  check_row(rates_death_other, which(
    rownames(rates_death_other$lambda_matrix) == "80-84 male never dead, other"
  ))
  check_row(rates_death_lung, which(
    rownames(rates_death_lung$lambda_matrix) == "distant 60-64 male 1"
  ))
})
