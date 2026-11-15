#' Step rates for death from other causes
#'
#' Piecewise-constant rates for the latent time to death from causes other than
#' lung cancer, in the form [nhppp::vdraw_sc_step()] expects. Derived from the
#' archival MILC cumulative-incidence tables `never.other`, `former.other` and
#' `current.other` by `data-raw/01_convert_mortality_other.R`.
#'
#' Each table \eqn{F} is the marginal CDF of an independent latent time (the
#' archival interpretation). On the table's knots \eqn{\tau_0 < \dots < \tau_J},
#' with \eqn{\bar F = 1 - F}, the rate on cell \eqn{j} is
#' \eqn{(\log \bar F(\tau_j) - \log \bar F(\tau_{j+1})) / (\tau_{j+1} - \tau_j)},
#' so the first arrival of the Poisson process reproduces \eqn{F} at every knot.
#' All strata are refined onto the union of their knots (a 0.25-year grid), which
#' leaves each step function unchanged. The clock is years since entry; the
#' stratum is fixed by the age group at entry for the whole follow-up.
#'
#' Rates cover the tabulated window only. With probability
#' `1 - F_max` there is no event in the window, and the generator applies the
#' tail rule. Where `terminal_certain` is `TRUE` (cumulative incidence reaches 1
#' at the last knot), the last cell's rate is set to 0 in `lambda_matrix` and the
#' generator draws a uniform time in `[terminal_lo, horizon]` for individuals with
#' no earlier event.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{time_breaks}{numeric vector of 65 breaks, 0 to 16 years since entry.}
#'   \item{lambda_matrix}{110 x 64 matrix of rates per year; row `i` is stratum
#'     `i` of `strata`; rownames are the archival table names.}
#'   \item{strata}{[data.table][data.table::data.table] with one row per stratum:
#'     `stratum_row`; `smoking_status` (`"never"`, `"former"`, `"current"`);
#'     `agegrp` (archival age-group label at entry, `"<40"`, ..., `">85"`);
#'     `sex`; `intensity` (`"1-10"`, `"11-20"`, `"21+"` for current smokers,
#'     `NA` otherwise); `horizon` (years); `F0` (cumulative incidence at 0);
#'     `F_max` (at the horizon); `n_knots`; `terminal_certain`; `terminal_lo`.}
#' }
#' @source Archival MILC 1.0 (Chrysanthopoulou, 2014), data `never.other`,
#'   `former.other`, `current.other`.
"rates_death_other"

#' Step rates for death from lung cancer
#'
#' Piecewise-constant rates for the latent time from diagnosis to death from
#' lung cancer, in the form [nhppp::vdraw_sc_step()] expects. Derived from the
#' archival MILC table `ci.lung` by
#' `data-raw/02_convert_mortality_lung_cancer.R`.
#'
#' Same construction as [rates_death_other]. The archival tables are in months;
#' breaks are converted to years since diagnosis and rates to per year. The
#' atom `F0` (cumulative incidence at diagnosis, up to 0.137 for distant stage) is
#' merged into the first cell, so the table is reproduced at every knot from the
#' first positive one on. The stratum is (stage at diagnosis, age group at
#' diagnosis, sex). Strata with age group `"NA"` are the archival fallback strata.
#' With probability `1 - F_max` there is no event in the window; the generator
#' then places lung-cancer death at `horizon` after diagnosis (archival behavior).
#'
#' @format A list with three elements:
#' \describe{
#'   \item{time_breaks}{numeric vector of 246 breaks, 0 to 251/12 years since
#'     diagnosis.}
#'   \item{lambda_matrix}{68 x 245 matrix of rates per year; row `i` is stratum
#'     `i` of `strata`; rownames are the archival table names.}
#'   \item{strata}{[data.table][data.table::data.table] with one row per stratum:
#'     `stratum_row`; `stage` (`"localized"`, `"regional"`, `"distant"`);
#'     `agegrp`; `sex`; `horizon` (years); `F0`; `F_max`; `n_knots`;
#'     `terminal_certain`; `terminal_lo`.}
#' }
#' @source Archival MILC 1.0 (Chrysanthopoulou, 2014), data `ci.lung`.
"rates_death_lung"
