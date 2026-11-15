# Helpers that convert the archival MILC cumulative-incidence (CIF) tables into
# the inputs nhppp::vdraw_sc_step() expects: one vector of time breaks and a
# matrix of piecewise-constant rates (one row per stratum).
#
# Model (model_description/model.tex, Section 2): each table F is the marginal
# CDF of an independent latent event time. On knots tau_0 < ... < tau_J with
# Fbar = 1 - F, the cell rate is
#   lambda_j = (log Fbar(tau_j) - log Fbar(tau_{j+1})) / (tau_{j+1} - tau_j),
# so the first arrival of the Poisson process reproduces F at every knot exactly.
# The lung atom F(0) > 0 is merged into the first cell by setting Fbar(tau_0) = 1.
#
# All strata of one table family are refined onto the union of their knots.
# A step function refined onto a finer grid is unchanged, so this is exact.

library(data.table)

#' Right-continuous CIF at its distinct knots.
#'
#' The archival tables list each step twice, (t, F-) and (t, F+); the value at
#' a knot is the largest estimate listed at that time.
#' @param cif list with elements `time`, `est` (archival format).
#' @return data.table with columns `time`, `F`, keyed by `time`.
cif_knots <- function(cif) {
  dt <- data.table(time = cif$time, est = cif$est)
  out <- dt[, .(F = max(est)), keyby = time]
  stopifnot(all(diff(out$F) >= 0), all(out$F >= 0), all(out$F <= 1))
  out
}

#' Piecewise-constant rates on a stratum's own knots.
#'
#' @param time knots tau_0 < ... < tau_J.
#' @param F CIF values at the knots.
#' @param merge_atom if TRUE, treat Fbar(tau_0) as 1, which merges the atom
#'   F(tau_0) into the first cell.
#' @return list: `lambda` (length J; the last cell may be Inf when F(tau_J) = 1),
#'   `terminal_certain` (TRUE when the last cell has infinite rate).
knot_rates <- function(time, F, merge_atom = TRUE) {
  J <- length(time) - 1L
  if (J < 1L) {
    return(list(lambda = numeric(0), terminal_certain = FALSE))
  }
  Fbar <- 1 - F
  if (merge_atom) {
    Fbar[1L] <- 1
  }
  lambda <- (log(Fbar[-(J + 1L)]) - log(Fbar[-1L])) / diff(time)
  # Cells after Fbar has reached 0 have no individuals at risk.
  lambda[is.nan(lambda)] <- 0
  inf_cells <- which(is.infinite(lambda))
  if (length(inf_cells) > 0L) {
    stopifnot(identical(inf_cells, J)) # only the last cell may be certain
  }
  list(lambda = lambda, terminal_certain = length(inf_cells) > 0L)
}

#' Refine a step function from its own knots onto a finer grid.
#'
#' @param time stratum knots (a subset of `grid`).
#' @param lambda rates on the stratum's cells (length(time) - 1).
#' @param grid union grid of breaks.
#' @return rates on the grid cells (length(grid) - 1); 0 beyond the stratum's
#'   last knot.
refine_to_grid <- function(time, lambda, grid) {
  out <- numeric(length(grid) - 1L)
  if (length(lambda) == 0L) {
    return(out)
  }
  stopifnot(all(time %in% grid))
  cell <- findInterval(grid[-length(grid)], time)
  inside <- cell >= 1L & cell <= length(lambda)
  out[inside] <- lambda[cell[inside]]
  out
}

#' Build the nhppp inputs for one family of tables.
#'
#' @param tables named list of archival CIF objects.
#' @param keys data.table with one row per element of `tables` (same order),
#'   holding the stratum key columns.
#' @param time_per_year table time units per year (1 for years, 12 for months).
#' @return list with
#'   - `time_breaks`: union grid of breaks, in years;
#'   - `lambda_matrix`: strata x cells matrix of rates per year. A certain
#'     terminal cell (F = 1 at the last knot) has rate 0 here and is flagged in
#'     `strata`;
#'   - `strata`: data.table with `stratum_row`, the key columns, `horizon`
#'     (years), `F0` (F at time 0), `F_max` (F at the last knot), `n_knots`,
#'     `terminal_certain`, `terminal_lo` (years; start of the certain cell).
build_rate_table <- function(tables, keys, time_per_year) {
  stopifnot(length(tables) == nrow(keys))
  knots <- lapply(tables, cif_knots)
  grid <- sort(unique(unlist(lapply(knots, `[[`, "time"))))

  rows <- vector("list", length(knots))
  lambda_matrix <- matrix(
    0,
    nrow = length(knots),
    ncol = length(grid) - 1L,
    dimnames = list(names(tables), NULL)
  )
  for (i in seq_along(knots)) {
    k <- knots[[i]]
    r <- knot_rates(k$time, k$F, merge_atom = TRUE)
    lam <- r$lambda
    terminal_lo <- NA_real_
    if (r$terminal_certain) {
      J <- length(lam)
      terminal_lo <- k$time[J] / time_per_year
      lam[J] <- 0
    }
    lambda_matrix[i, ] <- refine_to_grid(k$time, lam, grid) * time_per_year
    rows[[i]] <- data.table(
      horizon = max(k$time) / time_per_year,
      F0 = k$F[1L],
      F_max = max(k$F),
      n_knots = nrow(k),
      terminal_certain = r$terminal_certain,
      terminal_lo = terminal_lo
    )
  }
  strata <- cbind(
    data.table(stratum_row = seq_along(knots)),
    keys,
    rbindlist(rows)
  )
  list(
    time_breaks = grid / time_per_year,
    lambda_matrix = lambda_matrix,
    strata = strata[]
  )
}
