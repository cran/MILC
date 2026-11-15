# Tests for the canonical recipe synthetic_nlst_recipe and the recipe format.

test_that("synthetic_nlst_recipe has the recipe format and passes validation", {
  r <- synthetic_nlst_recipe
  expect_true(data.table::is.data.table(r))
  expect_identical(nrow(r), 26730L)
  expect_true(all(r$count == 1L))
  expect_equal(sum(r$weight), 1)
  expect_identical(
    sort(unique(r$smoking_status_at_spawn_age)), c("current", "former")
  )
  expect_true(all(r$study_entry_moment == as.Date("2002-06-30")))
  expect_identical(
    as.integer(format(r$inception_moment, "%Y")), r$birth_cohort
  )
  expect_no_error(validate_recipe(
    data.table::copy(r), drop_extra_variables = FALSE, drop_empty_strata = FALSE
  ))
})

test_that("synthetic smoking histories agree with status and NLST criteria", {
  r <- synthetic_nlst_recipe
  cur <- r[smoking_status_at_spawn_age == "current"]
  fmr <- r[smoking_status_at_spawn_age == "former"]
  expect_true(all(is.na(cur$smoking_quit_moment)))
  expect_true(all(cur$smoking_start_moment < cur$study_entry_moment))
  expect_true(all(fmr$smoking_start_moment < fmr$smoking_quit_moment))
  expect_true(all(fmr$smoking_quit_moment <= fmr$study_entry_moment))
  # NLST: former smokers quit within 15 years before entry
  yrs_quit <- as.numeric(fmr$study_entry_moment - fmr$smoking_quit_moment) / 365.25
  expect_true(all(yrs_quit <= 15))
  age_start <- as.numeric(r$smoking_start_moment - r$inception_moment) / 365.25
  expect_true(all(age_start >= 10 - 1 / 365))
})

test_that("validate_recipe rejects inconsistent smoking histories", {
  base <- synthetic_nlst_recipe[smoking_status_at_spawn_age == "former"][1]
  bad <- function(...) {
    x <- data.table::copy(base)
    x[, (...names()) := list(...)]
    x
  }
  expect_error(validate_recipe(bad(smoking_status_at_spawn_age = "smoker")), "must be")
  expect_error(validate_recipe(bad(smoking_quit_moment = as.Date(NA))), "former smokers")
  expect_error(
    validate_recipe(bad(smoking_quit_moment = base$study_entry_moment + 1)),
    "quit"
  )
  expect_error(
    validate_recipe(bad(inception_moment = base$study_entry_moment)),
    "birth_cohort"
  )
  expect_error(validate_recipe(bad(spawn_age = base$spawn_age + 3L)), "spawn_age")
  expect_error(
    validate_recipe(bad(
      smoking_status_at_spawn_age = "never", smoking_intensity_at_spawn_age = 0L
    )),
    "never smokers"
  )
  expect_error(validate_recipe(base[, !"smoking_start_moment"]), "missing")
})

test_that("every recipe row has an other-cause mortality table", {
  r <- synthetic_nlst_recipe
  rows <- death_other_stratum(
    r$spawn_age, r$sex, r$smoking_status_at_spawn_age,
    r$smoking_intensity_at_spawn_age
  )
  expect_false(anyNA(rows))
})

test_that("Population builds from the canonical recipe", {
  p <- suppressWarnings( # default rng_stream re-creates the RNGStreams singleton
    Population$new(1000, simulate_at_initialization = FALSE)
  )
  pop <- p$getPopulation()
  expect_identical(nrow(pop), 26730L)
  expect_true(all(c("study_entry_moment", "smoking_start_moment",
                    "smoking_quit_moment", "birth_cohort") %in% names(pop)))
  expect_identical(
    as.integer(format(pop$inception_moment, "%Y")), as.integer(pop$birth_cohort)
  )
})
