#' Validate population, (vectorized)
#'
#' @description Validate a population data.table. Meant for internal use.
#' Check that the population has the required variables (those of
#' [validate_recipe()] plus `person_id`, without `weight`) and, optionally,
#' that their values are consistent.
#'
#' Order, key, and indices
#'  This function orders the population by `person_id` and sets `person_id` to be the key.
#'  It also adds the index `birth_cohort__sex`
#'
#' @param population A data.table, data.frame, list, or matrix
#'               with the population data.  Coerced into a data.table
#' @param drop_extra_variables Logical. If TRUE, remove any variables that are
#' not in the required list
#' @param check Logical. If TRUE, perform checks.
#' @return A data.table with the validated `population`, ordered by the key `person_id`.
#'
#' @import data.table
validate_population <- function(
  population,
  drop_extra_variables = TRUE,
  check = TRUE
) {
  if (!is.data.table(population)) {
    population <- as.data.table(population)
  }
  pop_vars <- names(population)

  demographic_vars <- c(
    "person_id",
    "sex",
    "birth_cohort",
    "inception_moment",
    "study_entry_moment",
    "spawn_age"
  )
  risk_factor_vars <- c(
    "smoking_status_at_spawn_age",
    "smoking_intensity_at_spawn_age",
    "smoking_start_moment",
    "smoking_quit_moment"
  )

  required_vars <- c(demographic_vars, risk_factor_vars)

  missing_vars <- which(!required_vars %in% pop_vars)
  if (length(missing_vars) > 0) {
    stop(
      "These required variables are missing from `population`: ",
      paste(required_vars[missing_vars], collapse = ", "),
      ". See documentation for more information."
    )
  }

  if (drop_extra_variables) {
    drop_vars <- setdiff(pop_vars, required_vars)
    if (length(drop_vars)) population[, c(drop_vars) := NULL]
  }

  setkey(population, person_id)
  setorder(population, person_id)
  setindexv(population, c("birth_cohort", "sex"))

  if (isFALSE(check)) {
    return(population)
  }

  if (population[, anyNA(person_id) || anyDuplicated(person_id) > 0]) {
    stop("`person_id` must be non-missing and unique")
  }
  check_recipe_values(population)

  return(population)
}
