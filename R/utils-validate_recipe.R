#' Validate population recipe, (vectorized)
#'
#' @description Validate a population recipe for use in fastMILC
#' Check that the recipe file includes required variables,
#' and optionally perform checks
#'
#' Required variables:
#'  demographic_vars  <- c("sex", "birth_cohort")
#'  instruction_vars  <- c("weight", "spawn_age")
#'
#' Optional variables
#'  risk_factor_vars  <- c("smoking_status_at_NLST_entry", "smoking_intensity_at_NLST_entry")
#'
#' @param recipe A data.table, data.frame, list, or matrix
#'               with the population recipe. Coersed into a data.table
#' @param drop_extra_variables Logical. If TRUE, remove any variables that are
#' not in the required or optional lists above
#' @param drop_empty_strata Logical. If TRUE, remove any rows with missing or 0 `weight`
#' @param check Logical. If TRUE, perform checks.
#' @return A data.table with the validated `recipe`
#'
#' @import data.table

validate_recipe <- function(
  recipe,
  drop_extra_variables = TRUE,
  drop_empty_strata = TRUE,
  check = TRUE
) {
  if (!is.data.table(recipe)) {
    recipe <- as.data.table(recipe)
  }
  rx_vars <- names(recipe)

  # browser()

  demographic_vars <- c("sex", "birth_cohort")
  risk_factor_vars <- c(
    "smoking_status_at_spawn_age",
    "smoking_intensity_at_spawn_age"
  )
  instruction_vars <- c("weight", "spawn_age")

  required_vars <- c(demographic_vars, risk_factor_vars, instruction_vars)
  all_vars <- c(
    required_vars
  )

  missing_vars <- which(!required_vars %in% rx_vars)
  if (length(missing_vars) > 0) {
    stop(
      "These required variables are missing from `recipe`: ",
      paste(required_vars[missing_vars], collapse = ", "),
      ". See documentation for more information."
    )
  }

  if (drop_extra_variables) {
    drop_vars <- setdiff(rx_vars, all_vars)
    recipe[, c(drop_vars) := NULL]
  }
  if (drop_empty_strata) {
    recipe <- recipe[(!is.na(weight)) & weight > 0]
  }

  if (isFALSE(check)) {
    return(recipe)
  }

  stopifnot(recipe[, all(tolower(sex) %in% c("male", "female"))])
  stopifnot(recipe[, all(birth_cohort == floor(birth_cohort))])
  stopifnot(recipe[, all(
    spawn_age %between% list(0, getOption("milc.max_age"))
  )])

  return(recipe)
}
