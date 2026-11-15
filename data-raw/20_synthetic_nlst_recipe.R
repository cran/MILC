# Synthetic NLST-like population recipe: data-raw/nlst_recipe.rda -> package
# dataset `synthetic_nlst_recipe`, the canonical recipe for a simulation.
#
# The script puts the source into the recipe format that validate_recipe()
# expects:
#   - inception_moment = birth (mid-year of birth_cohort), the earliest moment
#     in the simulation: onset is simulated from birth;
#   - study_entry_moment = NLST entry (2002-06-30 in the source, where it was
#     stored as inception_moment); spawn_age = age at study entry;
#   - smoking status in the archival labels; every NLST participant who is not
#     a current smoker is a former smoker;
#   - smoking_start_moment, smoking_quit_moment: the source has no smoking
#     history, and the onset hazard switches at the ages of starting and
#     quitting (model.tex, Section 4), so the script draws them at random,
#     consistent with smoking status at entry and with the NLST inclusion
#     criterion that former smokers quit within 15 years before entry.
#
# The source has one row per stratum with a count. Smoking dates are drawn per
# person, so the recipe has one row per person (count = 1).
#
# Assumptions for the synthetic smoking history (not data):
#   - age at starting ~ Normal(START_MEAN, START_SD) truncated to
#     [START_MIN, age at quitting (former) or at entry (current) - MIN_YEARS];
#   - years from quitting to entry ~ Unif(0, QUIT_WINDOW) for former smokers.
#
# Run from the package root:
#   Rscript data-raw/20_synthetic_nlst_recipe.R

library(data.table)
pkgload::load_all(quiet = TRUE) # validate_recipe()

SEED <- 20021115L
START_MEAN <- 17
START_SD <- 4
START_MIN <- 10
MIN_YEARS <- 1 # minimum years smoked
QUIT_WINDOW <- 15 # NLST: former smokers quit within 15 years before entry
DAYS_PER_YEAR <- 365.25

env <- new.env()
load(file.path("data-raw", "nlst_recipe.rda"), envir = env)
src <- as.data.table(env$nlst_recipe)

stopifnot(
  setequal(
    names(src),
    c("sex", "birth_cohort", "spawn_age", "count", "weight",
      "smoking_status_at_spawn_age", "smoking_intensity_at_spawn_age",
      "inception_moment")
  ),
  !anyNA(src),
  src[, all(sex %in% c("male", "female"))],
  src[, all(count >= 1 & count == round(count))],
  src[, all(smoking_intensity_at_spawn_age >= 1)],
  # the source's inception_moment is the NLST entry date
  src[, all(spawn_age == as.integer(format(inception_moment, "%Y")) - birth_cohort)]
)

# One row per person ------------------------------------------------------
setorderv(
  src,
  c("sex", "smoking_status_at_spawn_age", "spawn_age", "birth_cohort",
    "smoking_intensity_at_spawn_age")
)
rcp <- src[rep(seq_len(.N), count)]

rcp <- rcp[, .(
  sex,
  birth_cohort = as.integer(birth_cohort),
  inception_moment = as.Date(sprintf("%d-06-30", as.integer(birth_cohort))),
  study_entry_moment = inception_moment,
  spawn_age = as.integer(spawn_age),
  # NLST enrolled current and former smokers only
  smoking_status_at_spawn_age = fifelse(
    smoking_status_at_spawn_age == "current_smoker", "current", "former"
  ),
  smoking_intensity_at_spawn_age = as.integer(smoking_intensity_at_spawn_age)
)]

# Synthetic smoking history ----------------------------------------------
set.seed(SEED)
n <- nrow(rcp)
former <- rcp$smoking_status_at_spawn_age == "former"

# age at quitting: within QUIT_WINDOW years before entry
age_quit <- rep(NA_real_, n)
age_quit[former] <- rcp$spawn_age[former] - QUIT_WINDOW * runif(sum(former))

# age at starting: truncated normal, ends at least MIN_YEARS before quitting
# (former) or entry (current)
upper <- fifelse(former, age_quit, rcp$spawn_age) - MIN_YEARS
stopifnot(all(upper > START_MIN))
p_lo <- pnorm(START_MIN, START_MEAN, START_SD)
p_hi <- pnorm(upper, START_MEAN, START_SD)
age_start <- qnorm(p_lo + runif(n) * (p_hi - p_lo), START_MEAN, START_SD)

rcp[, `:=`(
  smoking_start_moment = inception_moment + round(age_start * DAYS_PER_YEAR),
  smoking_quit_moment = inception_moment + round(age_quit * DAYS_PER_YEAR),
  count = 1L,
  weight = 1 / n
)]

synthetic_nlst_recipe <- rcp[]

# Checks against the recipe format ----------------------------------------
synthetic_nlst_recipe <- validate_recipe(
  synthetic_nlst_recipe,
  drop_extra_variables = FALSE,
  drop_empty_strata = FALSE,
  check = TRUE
)
stopifnot(
  identical(nrow(synthetic_nlst_recipe), as.integer(sum(src$count))),
  synthetic_nlst_recipe[former, all(
    smoking_quit_moment <= study_entry_moment &
      study_entry_moment - smoking_quit_moment <= round(QUIT_WINDOW * DAYS_PER_YEAR)
  )]
)

usethis::use_data(synthetic_nlst_recipe, overwrite = TRUE)
