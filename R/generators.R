# Vectorized event-time generators (nhppp-backed).
#
# SCAFFOLD STATUS: signatures and the intended nhppp calls are documented; bodies
# are stubs pending open decisions recorded in CLAUDE.md. Each returns a numeric
# vector of ages (or NA) aligned to the population rows.

#' Generate age at malignancy onset
#'
#' First arrival of the non-homogeneous Poisson process with cumulative intensity
#' [milc_HT_mal()], accounting for the piecewise smoking history
#' (never -> start age `ts` -> quit age `tq`). Intended implementation:
#' [nhppp::vdraw_cumulative_intensity()] with per-row `Lambda`/`Lambda_inv` closures
#' built from [milc_HT_mal()] (decision D2), `atmost1 = TRUE`.
#'
#' @param population A [data.table][data.table::data.table] with at least
#'   `sex`, smoking-history columns (`ts`, `tq`, intensity `d`), and the sampling
#'   window (`t_min`, `t_max`).
#' @return Numeric vector: age at onset, `NA` if none within the window.
#' @export
generate_onset <- function(population) {
  .not_yet_implemented("generate_onset", "D2 (cumulative-intensity onset sampler)")
}

#' Generate tumor progression and stage at diagnosis
#'
#' Gompertz tumor growth ([milc_t_vol()]) with random log-normal volume thresholds
#' at diagnosis / regional / distant, with rejection enforcing volume ordering and
#' bounds. Not a Poisson process; to be vectorized as a batched rejection sampler.
#'
#' @param population A [data.table][data.table::data.table] with onset ages and the
#'   log-normal volume parameters.
#' @return A [data.table][data.table::data.table] of `T_reg, T_dist, T_diagn, D_diagn, stage`.
#' @export
generate_progression <- function(population) {
  .not_yet_implemented("generate_progression", "porting t_prog() as a vectorized rejection sampler")
}

#' Generate age at death from other causes
#'
#' First arrival of the step-intensity process built from [other_cells] (the
#' converted marginal cumulative-incidence tables), stratified by entry
#' age-group x sex (x smoking intensity for current smokers). Intended
#' implementation: [nhppp::vdraw_sc_step_regular()] on a per-row rate matrix,
#' `atmost1 = TRUE`. Upper-tail handling (residual mass; F=1 strata) per the
#' pending decision in CLAUDE.md D1.
#'
#' @param population A [data.table][data.table::data.table] with `sex`, entry age,
#'   smoking `status`/intensity, and the sampling window.
#' @return Numeric vector: age at other-cause death.
#' @export
generate_death_other <- function(population) {
  .not_yet_implemented("generate_death_other", "D1 tail handling + rate-matrix assembly")
}

#' Generate age at death from lung cancer
#'
#' First arrival of the step-intensity process built from [lung_cells], stratified
#' by stage x age-group x sex. Same machinery as [generate_death_other()].
#'
#' @param population A [data.table][data.table::data.table] with `stage`, age at
#'   diagnosis, `sex`, and the sampling window.
#' @return Numeric vector: age at lung-cancer death (`NA`/censored per D1 tail rule).
#' @export
generate_death_cancer <- function(population) {
  .not_yet_implemented("generate_death_cancer", "D1 tail handling + rate-matrix assembly")
}

.not_yet_implemented <- function(fn, blocker) {
  stop(sprintf("%s() is a scaffold stub; blocked on: %s. See CLAUDE.md.", fn, blocker),
       call. = FALSE)
}
