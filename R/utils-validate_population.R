#' Validate population, (vectorized)
#'
#' @description Validate a population data.table, and create necessary placeholder variabels.
#' Meant for internal use.
#' Check that the population data.table includes required variables,
#' and optionally perform checks
#'
#' Required variables:
#'  demographic_vars <- c("person_id", "sex", "inception_moment", "spawn_age")
#'  risk_factor_vars <- c("smoking_status_at_spawn_age", "smoking_intensity_at_spawn_age")
#'
#' Order, key, and indices
#'  This function orders the population by `person_id` and sets `person_id` to be the key.
#'  It also adds the index `birth_cohort__sex`
#'
#' @param population A data.table, data.frame, list, or matrix
#'               with the population data.  Coersed into a data.table
#' @param drop_extra_variables Logical. If TRUE, remove any variables that are
#' not in the required or optional lists above
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
    "inception_moment",
    "spawn_age"
  )
  risk_factor_vars <- c(
    "smoking_status_at_spawn_age",
    "smoking_intensity_at_spawn_age"
  )

  required_vars <- c(demographic_vars, risk_factor_vars)
  all_vars <- c(required_vars)

  missing_vars <- which(!required_vars %in% pop_vars)
  if (length(missing_vars) > 0) {
    stop(
      "These required variables are missing from `population`: ",
      paste(required_vars[missing_vars], collapse = ", "),
      ". See documentation for more information."
    )
  }

  if (drop_extra_variables) {
    drop_vars <- setdiff(pop_vars, all_vars)
    population[, c(drop_vars) := NULL]
  }

  if (!"birth_cohort" %in% names(population)) {
    population[, birth_cohort := calendar_year_from_date(inception_moment)]
  } else {
    tmp_birth_cohort <- R.utils::tempvar(prefix = "birth_cohort_")
    population[,
      (tmp_birth_cohort) := calendar_year_from_date(inception_moment)
    ]
    stopifnot(population$birth_cohort == population[[tmp_birth_cohort]])
    population[, (tmp_birth_cohort) := NULL]
  }

  setkey(population, person_id)
  setorder(population, person_id)
  setindexv(population, c("birth_cohort", "sex"))

  if (isFALSE(check)) {
    return(population)
  }

  if (population[, any(is.na(person_id))]) {
    stop("Missing values in `person_id`")
  }
  stopifnot(population[, all(tolower(sex) %in% c("male", "female"))])
  stopifnot(population[, all(!is.na(inception_moment))])
  stopifnot(population[, all(
    spawn_age >= 0 & spawn_age <= getOption("milc.max_age")
  )])

  return(population)
}
