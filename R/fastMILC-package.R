#' fastMILC: Fast Vectorized MIcrosimulation Lung Cancer Model
#'
#' A vectorized re-implementation of the continuous-time microsimulation lung
#' cancer natural-history model MILC (Chrysanthopoulou, 2014). Event times are
#' drawn with the [nhppp](https://cran.r-project.org/package=nhppp) samplers
#' over a [data.table][data.table::data.table] population.
#'
#' @section Model components:
#' - Death from other causes ([generate_death_other()]) and from lung cancer
#'   ([generate_death_cancer()]): first arrival of a step-intensity process
#'   ([nhppp::vdraw_sc_step()]) built from the archival cumulative-incidence
#'   tables ([rates_death_other], [rates_death_lung]).
#' - Malignancy onset ([generate_onset()]): first arrival of the two-stage
#'   clonal expansion process over the smoking history, by
#'   [nhppp::vdraw_cumulative_intensity()] with a vectorized Newton inverse.
#' - Progression and stage ([generate_progression()]): Gompertz tumor growth
#'   with truncated log-normal volume thresholds.
#' - [Population]: an R6 class that simulates natural histories for a
#'   population expanded from a recipe ([synthetic_nlst_recipe]).
#'
#' The mathematics and the differences from the archival code are described in
#' `model_description/model.tex` in the source repository.
#'
#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @importFrom R6 R6Class
NULL
