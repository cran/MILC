# The shipped step-hazard cells encode the tabulated cumulative incidence
# CONDITIONAL on surviving the t=0 point mass F0 (the first-breakpoint CIF, which
# is > 0 for advanced-stage lung strata). Reconstruction must add the atom back:
#   F(t) = F0 + (1 - F0) * (1 - exp(-cumsum(lambda * width)))
# This must equal F at every breakpoint to fp precision.

roundtrip_ok <- function(cells, keys) {
  cells <- data.table::copy(cells)
  chk <- cells[, {
    F0 <- F_lo[1L]
    Lam <- cumsum(lambda * (t_hi - t_lo))
    F_recon <- F0 + (1 - F0) * (1 - exp(-Lam))
    .(max_abs_err = max(abs(F_recon - F_hi)))
  }, by = keys]
  chk[, max(max_abs_err)]
}

test_that("other_cells reproduce the tabulated CIF to fp precision", {
  err <- roundtrip_ok(other_cells, c("status", "agegrp", "sex", "intensity"))
  expect_lt(err, 1e-9)
})

test_that("lung_cells reproduce the tabulated CIF to fp precision", {
  err <- roundtrip_ok(lung_cells, c("stage", "agegrp", "sex"))
  expect_lt(err, 1e-9)
})

test_that("hazards are non-negative and cells are well-formed", {
  expect_true(other_cells[, all(t_hi > t_lo)])
  expect_true(lung_cells[,  all(t_hi > t_lo)])
  expect_true(other_cells[is.finite(lambda), all(lambda >= 0)])
  expect_true(lung_cells[is.finite(lambda),  all(lambda >= 0)])
})

test_that("infinite terminal hazards are confined to the documented F=1 strata", {
  # Only strata whose CIF reaches 1 may carry a non-finite (point-mass) hazard.
  bad_other <- other_cells[!is.finite(lambda), unique(.SD), .SDcols = c("status", "agegrp", "sex", "intensity")]
  expect_equal(nrow(bad_other), 1L)             # current >85 female 21+
  bad_lung <- lung_cells[!is.finite(lambda), unique(.SD), .SDcols = c("stage", "agegrp", "sex")]
  expect_equal(nrow(bad_lung), 1L)              # distant NA male (fallback)
})
