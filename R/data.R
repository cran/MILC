#' Other-cause death step-hazard cells
#'
#' Piecewise-constant step hazards for death from causes other than lung cancer,
#' derived from the archival MILC cumulative-incidence tables
#' (`never.other`/`former.other`/`current.other`) by the conversion in
#' `data-raw/01_convert_mortality_tables.R`. One row per (stratum, time cell) over
#' the tabulated horizon (time measured since entry into the stratum). The upper
#' tail (residual mass) is NOT included here; see CLAUDE.md decision D1.
#'
#' @format A [data.table][data.table::data.table] with columns:
#' \describe{
#'   \item{status}{smoking status: `"never"`, `"former"`, `"current"`.}
#'   \item{agegrp}{entry age group (`"<40"`, ..., `">85"`).}
#'   \item{sex}{`"female"` or `"male"`.}
#'   \item{intensity}{smoking intensity band for current smokers (`"1-10"`,
#'     `"11-20"`, `"21+"`), `NA` otherwise.}
#'   \item{t_lo, t_hi}{cell endpoints, years since entry.}
#'   \item{F_lo, F_hi}{cumulative incidence at the endpoints.}
#'   \item{lambda}{constant hazard on the cell (per year).}
#' }
#' @seealso [other_strata], [lung_cells]
"other_cells"

#' Other-cause death per-stratum metadata
#'
#' One row per stratum of [other_cells], carrying the information needed to append
#' upper-tail handling: tabulated horizon, cumulative-incidence plateau, residual
#' mass, cell count, and whether the stratum induces an infinite terminal hazard
#' (cumulative incidence reaching 1).
#'
#' @format A [data.table][data.table::data.table] with the stratum key columns plus
#' \describe{
#'   \item{F0}{point mass at t=0 (cumulative incidence at entry); 0 for all
#'     other-cause strata, `> 0` for advanced-stage lung strata. The step-hazard
#'     cells encode the distribution conditional on surviving this atom.}
#'   \item{horizon_yr}{tabulated horizon (years).}
#'   \item{F_max}{cumulative-incidence plateau.}
#'   \item{residual_mass}{`1 - F_max` (mass with no event in the tabulated window).}
#'   \item{n_cells}{number of step cells.}
#'   \item{any_inf_rate}{`TRUE` if the stratum's CIF reaches 1 (infinite terminal hazard).}
#' }
#' @seealso [other_cells]
"other_strata"

#' Lung-cancer death step-hazard cells
#'
#' Piecewise-constant step hazards for death from lung cancer, derived from the
#' archival `ci.lung` cumulative-incidence tables (time converted months to years).
#' One row per (stratum, time cell). Strata are keyed by stage x age-group x sex;
#' two fallback strata carry `agegrp = "NA"` (used when age group is unknown).
#'
#' @format A [data.table][data.table::data.table] with columns:
#' \describe{
#'   \item{stage}{`"localized"`, `"regional"`, `"distant"`.}
#'   \item{agegrp}{age group, or `"NA"` fallback.}
#'   \item{sex}{`"female"` or `"male"`.}
#'   \item{t_lo, t_hi}{cell endpoints, years since diagnosis.}
#'   \item{F_lo, F_hi}{cumulative incidence at the endpoints.}
#'   \item{lambda}{constant hazard on the cell (per year).}
#' }
#' @seealso [lung_strata], [other_cells]
"lung_cells"

#' Lung-cancer death per-stratum metadata
#'
#' One row per stratum of [lung_cells]; same fields as [other_strata].
#'
#' @format A [data.table][data.table::data.table] with the stratum key columns plus
#'   `horizon_yr`, `F_max`, `residual_mass`, `n_cells`, `any_inf_rate`.
#' @seealso [lung_cells]
"lung_strata"
