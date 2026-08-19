# data-raw/01_convert_mortality_tables.R
# Convert archival MILC CIF step-function tables (list-of-lists) into tidy
# data.tables of piecewise-constant step hazards over the TABULATED horizon.
#
# The upper-tail (residual-mass) handling is deliberately NOT applied here; it is
# decided separately. Each stratum's metadata (horizon, plateau F_max, residual
# mass) is retained so the tail can be appended later.
#
# Conversion (verified against tdeath_other/tdeath_lung inverse-CDF sampling):
#   latent survival  S(t) = 1 - F(t)          F = tabulated CIF (marginal/net)
#   step hazard      lambda_k = (log S(t_k) - log S(t_{k+1})) / (t_{k+1}-t_k)
# First arrival of the NHPPP with this step rate reproduces F at every breakpoint.

suppressMessages(library(data.table))

ARCHIVE <- "archival_code/data"

# ---- core converter: one CIF stratum -> (cells, summary) ---------------------
# el: list(time, est, var) with duplicated (t,F-)/(t,F+) step rows.
# time_scale: divide raw time by this (12 for ci.lung months -> years; 1 otherwise).
convert_stratum <- function(el, time_scale = 1) {
  browser()
  raw <- data.table(time = el$time / time_scale, est = el$est, var = el$var)
  raw[, interval := 2 - (1:.N) %% 2]
  tmp <- dcast(
    raw,
    est + var ~ interval,
    id.vars = c("est", "var"),
    value.var = "time"
  )
  setnames(tmp, old = names(tmp), new = c("F", "var_F", "t_a", "t_b"))
  tmp[, S := 1 - F]
  stopifnot(tmp[, all(diff(F) >= -1e-12)]) # F nondecreasing
  tmp[, dt := t_b - t_a]
  tmp[, lambda := -diff(log(S)) / dt] # >=0; Inf iff S hits 0
  # right-continuous CIF: post-jump value = max est at each unique time
  cdf <- raw[, .(F = max(est)), by = time][order(time)]
  cdf[, S := 1 - F]
  stopifnot(cdf[, all(diff(F) >= -1e-12)]) # F nondecreasing
  horizon <- cdf[, max(time)]
  F_max <- cdf[, max(F)]
  cells <- cdf[, .(
    t_lo = head(time, -1L),
    t_hi = time[-1L],
    F_lo = head(F, -1L),
    F_hi = F[-1L],
    lambda = -diff(log(S)) / diff(time) # >=0; Inf iff S hits 0
  )]
  summ <- data.table(
    F0 = cdf[1L, F], # point mass at t=0 (>0 for advanced-stage lung)
    horizon_yr = horizon,
    F_max = F_max,
    residual_mass = 1 - F_max,
    n_cells = cells[, .N],
    any_inf_rate = cells[, any(!is.finite(lambda))]
  )
  list(cells = cells, summary = summ)
}

# ---- apply over a named list of strata, parsing keys into covariate columns ----
# key_parser: function(key_string) -> named list of covariate columns.
convert_table <- function(obj, key_parser, time_scale = 1) {
  browser()
  parts <- lapply(names(obj), function(k) {
    covs <- as.data.table(key_parser(k))
    out <- convert_stratum(obj[[k]], time_scale = time_scale)
    # Replicate covs to the exact cell count so single-breakpoint strata (0 cells)
    # drop out instead of producing a recycled all-NA row.
    list(
      cells = cbind(covs[rep(1L, nrow(out$cells))], out$cells),
      summary = cbind(covs, out$summary)
    )
  })
  list(
    cells = rbindlist(lapply(parts, `[[`, "cells")),
    summary = rbindlist(lapply(parts, `[[`, "summary"))
  )
}

# ---- key parsers -------------------------------------------------------------
# other tables: "<agegrp> <sex> [<intensity> cigs] {never|former|current} dead, other"
parse_other <- function(key, status) {
  browser()
  tok <- strsplit(key, " ", fixed = TRUE)[[1]]
  if (status == "current") {
    list(status = "current", agegrp = tok[1], sex = tok[2], intensity = tok[3])
  } else {
    list(
      status = status,
      agegrp = tok[1],
      sex = tok[2],
      intensity = NA_character_
    )
  }
}
# ci.lung: "<stage> <agegrp> <sex> 1"   (agegrp may be literal "NA")
parse_lung <- function(key) {
  tok <- strsplit(key, " ", fixed = TRUE)[[1]]
  list(stage = tok[1], agegrp = tok[2], sex = tok[3])
}

# ---- run ---------------------------------------------------------------------
load(file.path(ARCHIVE, "never.other.rda"))
load(file.path(ARCHIVE, "former.other.rda"))
load(file.path(ARCHIVE, "current.other.rda"))
load(file.path(ARCHIVE, "ci.lung.rda"))

never <- convert_table(never.other, function(k) parse_other(k, "never"))
former <- convert_table(former.other, function(k) parse_other(k, "former"))
current <- convert_table(current.other, function(k) parse_other(k, "current"))
lung <- convert_table(ci.lung, parse_lung, time_scale = 12) # months -> years

# Combined other-cause death tables (status distinguishes the three sources)
other_cells <- rbindlist(list(never$cells, former$cells, current$cells))
other_strata <- rbindlist(list(never$summary, former$summary, current$summary))
lung_cells <- lung$cells
lung_strata <- lung$summary

setkey(other_cells, status, agegrp, sex, intensity, t_lo)
setkey(lung_cells, stage, agegrp, sex, t_lo)

# ---- persist (inspection artifacts; not yet package data/) -------------------
dir.create("data-raw/converted", showWarnings = FALSE, recursive = TRUE)
saveRDS(other_cells, "data-raw/converted/other_cells.rds")
saveRDS(other_strata, "data-raw/converted/other_strata.rds")
saveRDS(lung_cells, "data-raw/converted/lung_cells.rds")
saveRDS(lung_strata, "data-raw/converted/lung_strata.rds")

# ---- report ------------------------------------------------------------------
cat("== other_cells ==   rows:", nrow(other_cells), "\n")
print(head(other_cells, 4))
cat("\n== lung_cells ==    rows:", nrow(lung_cells), "\n")
print(head(lung_cells, 4))
cat("\n== other_strata (horizon / plateau / residual by status) ==\n")
print(other_strata[,
  .(
    n_strata = .N,
    horizon_yr = paste(sort(unique(horizon_yr)), collapse = ","),
    min_Fmax = round(min(F_max), 4),
    max_Fmax = round(max(F_max), 4),
    any_inf_rate = any(any_inf_rate)
  ),
  by = status
])
cat("\n== lung_strata (by stage) ==\n")
print(lung_strata[,
  .(
    n_strata = .N,
    horizon_yr = round(max(horizon_yr), 2),
    min_Fmax = round(min(F_max), 4),
    max_Fmax = round(max(F_max), 4),
    any_inf_rate = any(any_inf_rate)
  ),
  by = stage
])
cat("\nStrata whose CIF reaches F=1 (would give Inf terminal hazard):\n")
print(rbind(
  other_strata[F_max >= 1 - 1e-12],
  lung_strata[F_max >= 1 - 1e-12],
  fill = TRUE
))
