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

#' Synthetic NLST-like population recipe
#'
#' Synthetic NLST-like data: the canonical population recipe for a fastMILC
#' simulation, in the format [validate_recipe()] checks. One row per person:
#' 26,730 current and former smokers aged 49--79 who enter the study on
#' 2002-06-30 (the NLST enrollment period). `Population$new()` expands the
#' recipe to a population of the requested size in proportion to `weight`.
#' Built from `data-raw/nlst_recipe.rda` (one row per stratum with a count) by
#' `data-raw/20_synthetic_nlst_recipe.R`.
#'
#' The source has no smoking history. The onset hazard changes at the ages of
#' starting and quitting smoking, so the build script draws them at random
#' (fixed seed): age at starting from a Normal(17, 4) truncated to
#' \[10, age at quitting or at entry minus 1\]; years from quitting to entry
#' from Unif(0, 15) for former smokers (NLST inclusion criterion: quit within
#' 15 years). The NLST criterion of at least 30 pack-years is not imposed.
#'
#' @format A [data.table][data.table::data.table] with 26,730 rows and 11
#'   columns:
#' \describe{
#'   \item{sex}{`"male"` or `"female"`.}
#'   \item{birth_cohort}{integer, calendar year of birth (1923--1953).}
#'   \item{inception_moment}{`Date`, birth (June 30 of `birth_cohort`); the
#'     earliest moment in the simulation.}
#'   \item{study_entry_moment}{`Date`, study entry (2002-06-30).}
#'   \item{spawn_age}{integer, age at study entry (49--79 years).}
#'   \item{smoking_status_at_spawn_age}{`"current"` or `"former"`; NLST
#'     participants who are not current smokers are former smokers.}
#'   \item{smoking_intensity_at_spawn_age}{integer, cigarettes per day while
#'     smoking (10--201).}
#'   \item{smoking_start_moment}{`Date` of starting smoking (synthetic).}
#'   \item{smoking_quit_moment}{`Date` of quitting (synthetic); `NA` for
#'     current smokers.}
#'   \item{count}{integer, 1.}
#'   \item{weight}{`1 / 26730`.}
#' }
#' @source Synthetic data resembling the National Lung Screening Trial (NLST)
#'   population; `data-raw/nlst_recipe.rda`.
"synthetic_nlst_recipe"

#' Baseline model parameters
#'
#' One row per sex. The onset columns are the two-stage clonal expansion (TSCE)
#' parameters of the archival `HT_mal` (Chrysanthopoulou, 2014), which fastMILC
#' uses where they differ from Chrysanthopoulou (2017), Table 2
#' (`model_description/model.tex`, Table 1 and Section 9). Built by
#' `data-raw/10_baseline_model_parameters.R`.
#'
#' @format A [data.table][data.table::data.table] with 2 rows:
#' \describe{
#'   \item{sex}{`"male"`, `"female"`.}
#'   \item{n_stem_cells}{\eqn{X}, number of normal stem cells (\eqn{10^7}).}
#'   \item{normal_cell_initiation_rate}{\eqn{\nu_{ns}}.}
#'   \item{malignant_cell_conversion_rate}{\eqn{\mu}.}
#'   \item{cell_division_rate_nonsmoking}{\eqn{\alpha_{ns}}.}
#'   \item{cell_division_rate_smoking_base}{\eqn{\alpha_{s,0}}; while smoking
#'     \eqn{\alpha_s = \alpha_{s,0}\kappa(q)}.}
#'   \item{net_proliferation_rate}{\eqn{\gamma = \alpha - \beta - \mu} for
#'     non-smoking; while smoking \eqn{\gamma\kappa(q)}.}
#'   \item{smoking_initiation_multiplier}{\eqn{\rho}; \eqn{\nu_s = \rho\nu_{ns}}.}
#'   \item{smoking_dose_coefficient, smoking_dose_exponent}{\eqn{\alpha_1,
#'     \alpha_2} in \eqn{\kappa(q) = 1 + \alpha_1 q^{\alpha_2}}, \eqn{q}
#'     cigarettes per day.}
#'   \item{tumor_growth_rate}{\eqn{m}, Gompertz scale parameter (per day).}
#'   \item{volume_diagnosis_meanlog, volume_diagnosis_sdlog}{log-normal
#'     parameters of the tumor volume (ml) at diagnosis (archival `cdiagn`).}
#'   \item{volume_regional_meanlog, volume_regional_sdlog}{archival `creg`.}
#'   \item{volume_distant_meanlog, volume_distant_sdlog}{archival `cdist`. As
#'     in the archival code, [generate_progression()] draws the volume at the
#'     start of the regional stage with `cdist` and at the start of the
#'     distant stage with `creg`.}
#' }
#' Progression values are the posterior medians of Chrysanthopoulou (2017),
#' Table 2, with location = scale.
#' @source Archival MILC 1.0 (Chrysanthopoulou, 2014), `HT_mal`, `t_prog`.
"deterministic_params"
