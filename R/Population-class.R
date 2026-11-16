#' MILC microsimulation model (R6)
#'
#' Encapsulates a MILC natural-history simulation over a `data.table` population,
#' following the `nhppp`-DES idioms from the Trikalinos, Alarid,  Sereda, and Chrysanthopoulou
#' DES with `nhppp` paradigm.
#'
#' The simulation starts at study entry. Only people alive and without a
#' clinically diagnosed lung cancer at entry enter; undiagnosed malignancies at
#' entry are allowed. Simulated columns: `age_dead_from_other_causes` and
#' `age_dead_from_lung_cancer` (latent), `age_at_disease_onset`,
#' `age_at_regional_disease`, `age_at_distal_disease`, `age_at_diagnosis`,
#' `diameter_at_diagnosis`, `stage_at_diagnosis` (recorded only before the
#' end of follow-up), `age_dead`, `cause_of_death` (`NA` if alive at
#' `getOption("milc.max_age")`).
#'
#' @export
Population <- R6::R6Class(
  classname = "Population",
  public = list(
    #' @description Build the population from a recipe and, optionally,
    #'   simulate it.
    #' @param target_population_size number of people, a positive integer.
    #'   People are drawn with replacement from the recipe rows with
    #'   probability proportional to `weight` (multinomial counts per row, from
    #'   the `"population"` stream), so the population has exactly this size.
    #' @param recipe population recipe; see [validate_recipe()]. Default
    #'   [synthetic_nlst_recipe].
    #' @param parameters model parameters, one row per sex; default
    #'   [deterministic_params].
    #' @param rng_stream an [RNGStreams] object (one stream per model
    #'   component).
    #' @param simulate_at_initialization if `TRUE`, run `update()`.
    initialize = function(
      target_population_size,
      recipe = fastMILC::synthetic_nlst_recipe,
      parameters = fastMILC::deterministic_params,
      rng_stream = RNGStreams$new(),
      simulate_at_initialization = TRUE
    ) {
      private$.rng_stream <- rng_stream
      private$.rng_stream$resetAllStreams()
      private$.rng_stream$activateStream("population")
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )

      if (!is.numeric(target_population_size) ||
          length(target_population_size) != 1L ||
          is.na(target_population_size) ||
          target_population_size < 1 ||
          target_population_size != round(target_population_size)) {
        stop("`target_population_size` must be a positive integer")
      }
      private$.target_population_size <- as.integer(target_population_size)
      private$.parameters <- parameters
      private$.generatePopulationFromRecipe(
        recipe = recipe,
        drop_extra_variables = FALSE,
        drop_empty_strata = TRUE,
        check = TRUE
      )
      private$.updatePersonParameterValues(reset_rng_stream = TRUE)
      private$.simulation_erased <- TRUE

      if (simulate_at_initialization) {
        self$update()
      }
      invisible(self)
    },

    #' @description Run the natural-history simulation, filling event columns.
    update = function() {
      stopifnot(private$.simulation_erased)
      private$.simulated_cols <- NULL
      private$.generateAgeDeadFromOtherCauses(reset_rng_stream = TRUE)
      private$.generateDiseaseOnset(reset_rng_stream = TRUE)
      private$.generateDiseaseProgression(reset_rng_stream = TRUE)
      private$.enforceEntryCondition()
      private$.generateAgeDeadFromLungCancer(reset_rng_stream = TRUE)
      private$.assembleHistory()
      private$.simulation_erased <- FALSE
      invisible(self)
    },

    #' @description Clear simulated columns so the model can be re-run.
    reset = function() {
      drop_cols <- intersect(
        names(private$.population),
        private$.simulated_cols
      )
      if (length(drop_cols)) {
        private$.population[, (drop_cols) := NULL]
      }
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
    .target_population_size = NULL,
    .population_size = NULL,
    .recipe = NULL,
    .simulation_erased = NULL,
    .simulated_cols = NULL,
    .rng_stream = NULL,
    .parameters = NULL,

    # ' @description return a validated recipe object
    # ' @param recipe A data.table, data.frame, list, or matrix
    # '               with the population recipe. Coersed into a data.table
    # ' @param drop_extra_variables Logical. If TRUE, remove any variables that are
    # '                             not in the required or optional lists above
    # ' @param drop_empty_strata Logical. If TRUE, remove any rows with missing or 0 `weight`
    # ' @param check Logical. If TRUE, perform checks.
    .validateRecipe = function(
      recipe,
      drop_extra_variables,
      drop_empty_strata,
      check
    ) {
      return(
        validate_recipe(
          recipe = recipe,
          drop_extra_variables = drop_extra_variables,
          drop_empty_strata = drop_empty_strata,
          check = check
        )
      )
    },

    # ' @description Parse and check the recipe file. See `fastMILC::validateRecipe()`
    # ' @param recipe A data.table, data.frame, list, or matrix
    # '               with the population recipe. Coersed into a data.table
    # ' @param drop_extra_variables Logical. If TRUE, remove any variables in the recipe that are
    # '                             not in the required or optional lists above
    # ' @param drop_empty_strata Logical. If TRUE, remove any rows in the recipe with missing or 0 `weight`
    # ' @param check Logical. If TRUE, perform checks.
    .generatePopulationFromRecipe = function(
      recipe,
      drop_extra_variables,
      drop_empty_strata,
      check
    ) {
      if (!is.data.table(recipe)) {
        setDT(recipe)
      }
      recipe <- copy(recipe)
      private$.recipe <- private$.validateRecipe(
        recipe = recipe,
        drop_extra_variables = drop_extra_variables,
        drop_empty_strata = drop_empty_strata,
        check = check
      )

      # probability sampling: multinomial counts per recipe row, total exactly
      # target_population_size (the "population" stream is active)
      private$.recipe[, `:=`(
        sim_count = as.vector(stats::rmultinom(
          1L, private$.target_population_size, weight
        )),
        stratum_id = seq_len(.N)
      )]

      recipe_expanded <- private$.recipe[rep(stratum_id, sim_count)]
      recipe_expanded[, `:=`(
        person_id = seq_len(.N)
      )]

      # .validatePopulation sets key `person_id` and indices
      private$.population <-
        private$.validatePopulation(
          population = recipe_expanded,
          drop_extra_variables = TRUE,
          check = TRUE
        )

      private$.population_size <- nrow(private$.population)
      invisible(self)
    },

    .validatePopulation = function(
      population = recipe_expanded,
      drop_extra_variables = TRUE,
      check = TRUE
    ) {
      validate_population(
        population = population,
        drop_extra_variables = drop_extra_variables,
        check = check
      )
    },

    .updatePersonParameterValues = function(reset_rng_stream) {
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )
      if (reset_rng_stream) {
        private$.rng_stream$resetStream("parameters")
      }
      private$.rng_stream$activateStream("parameters")

      # One parameter row per sex (data-raw/10_baseline_model_parameters.R).
      private$.population <- private$.parameters[
        private$.population,
        on = "sex"
      ]
      setkey(private$.population, person_id)

      # Ages at starting and quitting smoking (NA if never / not quit).
      private$.population[, `:=`(
        age_smoking_start = as.numeric(
          smoking_start_moment - inception_moment
        ) / 365.25,
        age_smoking_quit = as.numeric(
          smoking_quit_moment - inception_moment
        ) / 365.25
      )]

      invisible(self)
    },

    # Columns of the parameter table used by the generators (all but `sex`).
    .parameterCols = function() {
      setdiff(names(private$.parameters), "sex")
    },

    # Latent age at death from other causes (model.tex, Section 6). Draws with
    # death before study entry (possible only under the archival tail for
    # entry > milc.max_age - horizon) are redrawn: only people alive at entry
    # enter the simulation.
    .generateAgeDeadFromOtherCauses = function(reset_rng_stream,
                                               max_rounds = 1000L) {
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )
      if (reset_rng_stream) {
        private$.rng_stream$resetStream("death_from_other_causes")
      }
      private$.rng_stream$activateStream("death_from_other_causes")

      draw <- function(i) {
        private$.population[i, generate_death_other(
          age_entry = spawn_age,
          sex = sex,
          smoking_status = smoking_status_at_spawn_age,
          intensity = smoking_intensity_at_spawn_age,
          tail = getOption("milc.death_other_tail")
        )]
      }
      vals <- draw(seq_len(nrow(private$.population)))
      private$.population[, age_dead_from_other_causes := vals]
      for (round in seq_len(max_rounds)) {
        bad <- private$.population[, which(age_dead_from_other_causes < spawn_age)]
        if (length(bad) == 0L) break
        set(private$.population, bad, "age_dead_from_other_causes", draw(bad))
      }
      if (length(bad) > 0L) {
        stop("Other-cause death before entry not resolved after ",
             max_rounds, " rounds for ", length(bad), " individual(s)")
      }

      private$.simulated_cols <- c(
        private$.simulated_cols,
        "age_dead_from_other_causes"
      )
      invisible(self)
    },

    # Latent age at onset of the first malignant cell, from birth
    # (model.tex, Section 4). NA if no onset by milc.max_age.
    .drawOnset = function(rows) {
      private$.population[
        rows,
        generate_onset(
          params = .SD,
          age_smoking_start = age_smoking_start,
          age_smoking_quit = age_smoking_quit,
          intensity = smoking_intensity_at_spawn_age
        ),
        .SDcols = private$.parameterCols()
      ]
    },

    .generateDiseaseOnset = function(reset_rng_stream) {
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )
      if (reset_rng_stream) {
        private$.rng_stream$resetStream("disease_onset")
      }
      private$.rng_stream$activateStream("disease_onset")

      vals <- private$.drawOnset(seq_len(nrow(private$.population)))
      private$.population[, age_at_disease_onset := vals]

      private$.simulated_cols <- c(
        private$.simulated_cols,
        "age_at_disease_onset"
      )
      invisible(self)
    },

    # Latent progression after onset (model.tex, Section 5).
    .progressionCols = c(
      "age_at_regional_disease",
      "age_at_distal_disease",
      "age_at_diagnosis",
      "diameter_at_diagnosis",
      "stage_at_diagnosis"
    ),

    .drawProgression = function(rows) {
      generate_progression(
        age_onset = private$.population$age_at_disease_onset[rows],
        params = private$.population[
          rows, .SD,
          .SDcols = private$.parameterCols()
        ]
      )
    },

    .generateDiseaseProgression = function(reset_rng_stream) {
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )
      if (reset_rng_stream) {
        private$.rng_stream$resetStream("disease_progression")
      }
      private$.rng_stream$activateStream("disease_progression")

      prog <- private$.drawProgression(seq_len(nrow(private$.population)))
      private$.population[, (private$.progressionCols) := prog]

      private$.simulated_cols <- c(
        private$.simulated_cols,
        private$.progressionCols
      )
      invisible(self)
    },

    # Entry condition: the simulation starts at study entry (trial start), and
    # only people without a clinically diagnosed lung cancer at entry enter.
    # Undiagnosed (pre-clinical) malignancies at entry are allowed. Rejection
    # sampling: for people diagnosed before entry, onset and progression are
    # redrawn jointly until diagnosis is at or after entry (or does not occur),
    # which samples their conditional law given eligibility exactly. Not in
    # the archival code (model.tex, Section 2).
    .enforceEntryCondition = function(max_rounds = 1000L) {
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )
      pop <- private$.population # set() on existing columns only
      for (round in seq_len(max_rounds)) {
        bad <- pop[, which(!is.na(age_at_diagnosis) & age_at_diagnosis < spawn_age)]
        if (length(bad) == 0L) break
        private$.rng_stream$activateStream("disease_onset")
        set(pop, bad, "age_at_disease_onset", private$.drawOnset(bad))
        private$.rng_stream$activateStream("disease_progression")
        prog <- private$.drawProgression(bad)
        for (v in private$.progressionCols) set(pop, bad, v, prog[[v]])
      }
      if (length(bad) > 0L) {
        stop("Entry condition not met after ", max_rounds, " rounds for ",
             length(bad), " individual(s)")
      }
      invisible(self)
    },

    # Latent age at death from lung cancer, for everyone with a (latent)
    # diagnosis, as in the archival nat_hist (model.tex, Section 7).
    .generateAgeDeadFromLungCancer = function(reset_rng_stream) {
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )
      if (reset_rng_stream) {
        private$.rng_stream$resetStream("death_from_lung_cancer")
      }
      private$.rng_stream$activateStream("death_from_lung_cancer")

      private$.population[,
        age_dead_from_lung_cancer := generate_death_cancer(
          age_diagnosis = age_at_diagnosis,
          stage = stage_at_diagnosis,
          sex = sex
        )
      ]

      private$.simulated_cols <- c(
        private$.simulated_cols,
        "age_dead_from_lung_cancer"
      )
      invisible(self)
    },

    # Assembly (model.tex, Section 8). Death is the first of the two latent
    # death times (cause "other_causes" if the other-cause time is earlier,
    # "lung_cancer" otherwise, as in nat_hist). Follow-up ends at death or at
    # milc.max_age; a person alive at milc.max_age has age_dead and
    # cause_of_death NA. Progression events (onset, regional, distant,
    # diagnosis) after the end of follow-up are set to NA, together with stage
    # and diameter at diagnosis; onset and stage transitions before entry are
    # kept (undiagnosed malignancy at entry). The cause-specific death ages
    # stay latent.
    .assembleHistory = function() {
      M <- getOption("milc.max_age")
      private$.population[, age_dead := pmin(
        age_dead_from_other_causes, age_dead_from_lung_cancer,
        na.rm = TRUE
      )]
      private$.population[, cause_of_death := fifelse(
        !is.na(age_dead_from_lung_cancer) &
          age_dead_from_lung_cancer <= age_dead_from_other_causes,
        "lung_cancer", "other_causes"
      )]
      pop <- private$.population # set() on existing columns only
      pop[age_dead > M, `:=`(age_dead = NA_real_, cause_of_death = NA_character_)]
      end_of_followup <- pop[, fifelse(is.na(age_dead), M, age_dead)]
      for (v in c("age_at_disease_onset", "age_at_regional_disease",
                  "age_at_distal_disease", "age_at_diagnosis")) {
        after <- which(pop[[v]] > end_of_followup)
        if (length(after)) set(pop, after, v, NA_real_)
      }
      pop[is.na(age_at_diagnosis), `:=`(
        diameter_at_diagnosis = NA_real_,
        stage_at_diagnosis = NA_character_
      )]
      private$.simulated_cols <- c(
        private$.simulated_cols,
        "age_dead",
        "cause_of_death"
      )
      invisible(self)
    }
  )
)
