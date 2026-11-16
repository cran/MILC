# End-to-end tests of the Population class on the canonical recipe.

get_sim <- local({
  sim <- NULL
  function() {
    if (is.null(sim)) {
      # the default rng_stream re-creates the RNGStreams singleton (warns)
      p <- suppressWarnings(Population$new(1000))
      sim <<- list(p = p, pop = p$getPopulation())
    }
    sim
  }
})

test_that("only people alive and undiagnosed at entry enter", {
  pop <- get_sim()$pop
  expect_false(any(pop$age_at_diagnosis < pop$spawn_age, na.rm = TRUE))
  expect_false(any(pop$age_dead < pop$spawn_age, na.rm = TRUE))
  expect_false(any(pop$age_dead_from_other_causes < pop$spawn_age))
  # undiagnosed malignancies at entry are allowed and occur
  expect_true(any(pop$age_at_disease_onset < pop$spawn_age, na.rm = TRUE))
})

test_that("assembled histories are internally consistent", {
  pop <- get_sim()$pop
  M <- getOption("milc.max_age")
  with(pop, {
    end <- ifelse(is.na(age_dead), M, age_dead)
    expect_true(all(age_dead == pmin(age_dead_from_other_causes,
                                     age_dead_from_lung_cancer, na.rm = TRUE),
                    na.rm = TRUE))
    lung <- cause_of_death == "lung_cancer" & !is.na(cause_of_death)
    expect_true(all(age_dead[lung] == age_dead_from_lung_cancer[lung]))
    expect_true(all(!is.na(age_at_diagnosis[lung])))
    for (v in list(age_at_disease_onset, age_at_regional_disease,
                   age_at_distal_disease, age_at_diagnosis)) {
      expect_true(all(v <= end, na.rm = TRUE))
    }
    expect_identical(is.na(stage_at_diagnosis), is.na(age_at_diagnosis))
    expect_identical(is.na(diameter_at_diagnosis), is.na(age_at_diagnosis))
    dx <- !is.na(age_at_diagnosis)
    expect_true(all(age_at_diagnosis[dx] >= age_at_disease_onset[dx]))
    expect_true(all(age_dead_from_lung_cancer[dx] >= age_at_diagnosis[dx]))
    expect_true(all(!is.na(age_at_disease_onset[lung])))
  })
})

test_that("reset() and update() reproduce the simulation", {
  sim <- get_sim()
  sim$p$reset()
  sim$p$update()
  expect_identical(sim$p$getPopulation(), sim$pop)
})

test_that("probability sampling gives the target size, follows weights, and is reproducible", {
  expect_identical(nrow(get_sim()$pop), 1000L)
  rec <- data.table::copy(synthetic_nlst_recipe)
  cur <- rec$smoking_status_at_spawn_age == "current"
  rec[, weight := ifelse(cur, 3, 1)]
  target <- 3 * sum(cur) / (3 * sum(cur) + sum(!cur))
  n <- 20000L
  draw <- function() {
    suppressWarnings(Population$new(n, recipe = rec,
                                    simulate_at_initialization = FALSE))$getPopulation()
  }
  pop1 <- draw()
  expect_identical(nrow(pop1), n)
  share <- mean(pop1$smoking_status_at_spawn_age == "current")
  expect_lt(abs(share - target), 5 * sqrt(target * (1 - target) / n))
  expect_identical(draw(), pop1)
  # a target above the recipe size samples rows more than once
  expect_identical(nrow(suppressWarnings(
    Population$new(30000L, simulate_at_initialization = FALSE)
  )$getPopulation()), 30000L)
})

test_that("target_population_size is validated", {
  for (bad in list(0, -1, 10.5, NA, c(10, 20), "10")) {
    expect_error(suppressWarnings(
      Population$new(bad, simulate_at_initialization = FALSE)
    ), "positive integer")
  }
})
