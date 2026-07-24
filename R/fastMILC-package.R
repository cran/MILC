#' fastMILC: Fast Vectorized MIcrosimulation Lung Cancer Model
#'
#' A vectorized re-implementation of the continuous-time microsimulation lung
#' cancer natural-history model MILC. Event times are drawn with the
#' [nhppp](https://cran.r-project.org/package=nhppp) samplers over a
#' [data.table][data.table::data.table] population.
#'
#' @section Model components:
#' - Malignancy onset (`generate_onset`): first arrival of a non-homogeneous
#'   Poisson process with the closed-form two-stage clonal-expansion cumulative
#'   intensity [milc_HT_mal()], drawn via [nhppp::vdraw_cumulative_intensity()].
#' - Tumor progression and stage (`generate_progression`): Gompertz tumor growth
#'   ([milc_t_vol()]) with random log-normal volume thresholds.
#' - Death from other causes (`generate_death_other`) and from lung cancer
#'   (`generate_death_cancer`): first arrival of a step-intensity process built
#'   from the tabulated cumulative-incidence data ([other_cells], [lung_cells]).
#'
#' @keywords internal
"_PACKAGE"

## data.table awareness: suppress R CMD check NOTEs for NSE column symbols.
## (Kept minimal; extend as generators are implemented.)
utils::globalVariables(c(
  ".", "time", "est", "F", "S", "F_lo", "F_hi", "t_lo", "t_hi", "lambda",
  "status", "agegrp", "sex", "intensity", "stage"
))

#' @importFrom R6 R6Class
#' @import data.table
NULL
