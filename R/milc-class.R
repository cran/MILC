#' MILC microsimulation model (R6)
#'
#' Encapsulates a MILC natural-history simulation over a `data.table` population,
#' following the `nhppp`-DES idioms from the Trikalinos, Alarid,  Sereda, and Chrysanthopoulou  
#' DES with `nhppp` paradigm. 
#' object, event generators fill columns in place, and results are retrieved with
#' [getPopulation()][MILC]. Sampling uses the `nhppp` package.
#'
#' @export
Population <- R6::R6Class(
  classname = "Population",
  public = list(
    initialize = function(
      
    ) {
      private$.setPopulation(population = population, copy = TRUE)
      private$.simulation_erased <- TRUE
      invisible(self)
    },

    #' @description Run the natural-history simulation, filling event columns.
    update = function() {
      stopifnot(private$.simulation_erased)
      private$.generate_death_other()
      private$.generate_onset()
      private$.generate_progression()
      private$.generate_death_cancer()
      private$.resolve_trajectory()
      private$.simulation_erased <- FALSE
      invisible(self)
    },

    #' @description Clear simulated columns so the model can be re-run.
    reset = function() {
      drop_cols <- intersect(names(private$.population), private$.simulated_cols)
      if (length(drop_cols)) private$.population[, (drop_cols) := NULL]
      private$.simulation_erased <- TRUE
      invisible(self)
    },

    #' @description Return the population table.
    #' @param copy If `TRUE` (default) return a deep copy; if `FALSE` a reference.
    getPopulation = function(copy = TRUE) {
      if (copy) copy(private$.population) else private$.population
    }
  ),
  private = list(
    .population = NULL,
    .simulation_erased = NULL,
    .simulated_cols = c(
      "age_onset", "T_reg", "T_dist", "T_diagn", "D_diagn", "stage",
      "age_death_other", "age_death_cancer", "age_death", "cause", "with_cancer"
    ),

    .setPopulation = function(population, copy) {
      stopifnot(data.table::is.data.table(population))
      required <- c("id", "sex", "spawn_age", "max_simulation_age")
      stopifnot(all(required %in% names(population)))
      private$.population <- if (copy) data.table::copy(population) else population
      data.table::setindex(private$.population, id)
      invisible(self)
    },

    .generate_death_other  = function() generate_death_other(private$.population),
    .generate_onset        = function() generate_onset(private$.population),
    .generate_progression  = function() generate_progression(private$.population),
    .generate_death_cancer = function() generate_death_cancer(private$.population),

    .resolve_trajectory = function() {
      .not_yet_implemented(".resolve_trajectory",
                           "stitching min(other, cancer) + incidence/exclusion logic (nat_hist)")
    }
  )
)
