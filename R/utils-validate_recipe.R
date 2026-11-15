#' Validate population recipe, (vectorized)
#'
#' @description Validate a population recipe for use in fastMILC. Check that
#' the recipe has the required variables and, optionally, that their values are
#' consistent. [synthetic_nlst_recipe] is the canonical recipe.
#'
#' Required variables:
#' \describe{
#'   \item{`sex`}{`"male"` or `"female"`.}
#'   \item{`birth_cohort`}{integer calendar year of birth.}
#'   \item{`inception_moment`}{`Date`, birth: the earliest moment in the
#'     simulation (onset is simulated from birth). Its calendar year is
#'     `birth_cohort`.}
#'   \item{`study_entry_moment`}{`Date`, entry into the study.}
#'   \item{`spawn_age`}{age at study entry (years); within one year of
#'     `study_entry_moment - inception_moment`.}
#'   \item{`smoking_status_at_spawn_age`}{`"never"`, `"former"`, or
#'     `"current"`.}
#'   \item{`smoking_intensity_at_spawn_age`}{cigarettes per day while smoking;
#'     0 for never smokers, at least 1 otherwise.}
#'   \item{`smoking_start_moment`}{`Date` of starting smoking; `NA` for never
#'     smokers; between `inception_moment` and `study_entry_moment` otherwise.}
#'   \item{`smoking_quit_moment`}{`Date` of quitting; non-missing for former
#'     smokers only, after `smoking_start_moment` and not after
#'     `study_entry_moment`.}
#'   \item{`weight`}{non-negative sampling weight.}
#' }
#' Optional variables: `count`.
#'
#' @param recipe A data.table, data.frame, list, or matrix
#'               with the population recipe. Coerced into a data.table
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

  demographic_vars <- c(
    "sex", "birth_cohort", "inception_moment", "study_entry_moment"
  )
  risk_factor_vars <- c(
    "smoking_status_at_spawn_age",
    "smoking_intensity_at_spawn_age",
    "smoking_start_moment",
    "smoking_quit_moment"
  )
  instruction_vars <- c("weight", "spawn_age")
  optional_vars <- c("count")

  required_vars <- c(demographic_vars, risk_factor_vars, instruction_vars)
  all_vars <- c(required_vars, optional_vars)

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
    if (length(drop_vars)) recipe[, c(drop_vars) := NULL]
  }
  if (drop_empty_strata) {
    recipe <- recipe[(!is.na(weight)) & weight > 0]
  }

  if (isFALSE(check)) {
    return(recipe)
  }

  check_recipe_values(recipe)
  return(recipe)
}

#' Value checks shared by validate_recipe() and validate_population()
#'
#' @param dt data.table with the recipe variables.
#' @return `invisible(TRUE)`; stops with a message naming the failed check.
#' @keywords internal
check_recipe_values <- function(dt) {
  fail <- function(msg, rows) {
    stop(msg, " (", sum(rows), " row(s); first: ", which(rows)[1L], ")",
         call. = FALSE)
  }
  chk <- function(ok, msg) {
    ok[is.na(ok)] <- FALSE
    if (!all(ok)) fail(msg, !ok)
  }
  for (v in c("inception_moment", "study_entry_moment",
              "smoking_start_moment", "smoking_quit_moment")) {
    if (!inherits(dt[[v]], "Date")) stop("`", v, "` must be a Date", call. = FALSE)
  }
  chk(dt$sex %in% c("male", "female"), "`sex` must be \"male\" or \"female\"")
  chk(dt$birth_cohort == floor(dt$birth_cohort), "`birth_cohort` must be integer")
  chk(
    calendar_year_from_date(dt$inception_moment) == dt$birth_cohort,
    "the year of `inception_moment` (birth) must equal `birth_cohort`"
  )
  chk(
    dt$study_entry_moment >= dt$inception_moment,
    "`study_entry_moment` must not precede `inception_moment`"
  )
  chk(
    dt$spawn_age >= 0 & dt$spawn_age <= getOption("milc.max_age"),
    "`spawn_age` must be in [0, milc.max_age]"
  )
  age_at_entry <- as.numeric(dt$study_entry_moment - dt$inception_moment) / 365.25
  chk(
    abs(dt$spawn_age - age_at_entry) < 1,
    "`spawn_age` must be the age at `study_entry_moment`"
  )

  status <- dt$smoking_status_at_spawn_age
  chk(
    status %in% c("never", "former", "current"),
    "`smoking_status_at_spawn_age` must be \"never\", \"former\" or \"current\""
  )
  q <- dt$smoking_intensity_at_spawn_age
  chk(
    fifelse(status == "never", q == 0, q >= 1),
    "`smoking_intensity_at_spawn_age` must be 0 for never smokers and >= 1 otherwise"
  )
  start <- dt$smoking_start_moment
  quit <- dt$smoking_quit_moment
  never <- status == "never"
  former <- status == "former"
  current <- status == "current"
  chk(!never | (is.na(start) & is.na(quit)),
      "never smokers must have no smoking start or quit moment")
  chk(!current | (!is.na(start) & is.na(quit)),
      "current smokers must have a start moment and no quit moment")
  chk(!former | (!is.na(start) & !is.na(quit)),
      "former smokers must have start and quit moments")
  smoker <- !never
  chk(
    !smoker | (start >= dt$inception_moment & start <= dt$study_entry_moment),
    "`smoking_start_moment` must be between birth and study entry"
  )
  chk(
    !former | (quit > start & quit <= dt$study_entry_moment),
    "`smoking_quit_moment` must be after starting and not after study entry"
  )
  if ("weight" %in% names(dt)) {
    chk(is.na(dt$weight) | dt$weight >= 0, "`weight` must be non-negative")
  }
  invisible(TRUE)
}
