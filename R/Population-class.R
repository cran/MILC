#' MILC microsimulation model (R6)
#'
#' Encapsulates a MILC natural-history simulation over a `data.table` population,
#' following the `nhppp`-DES idioms from the Trikalinos, Alarid,  Sereda, and Chrysanthopoulou
#' DES with `nhppp` paradigm.
#'
#' @export
Population <- R6::R6Class(
  classname = "Population",
  public = list(
    initialize = function(
      target_population_size,
      recipe,
      parameters = fastMILC::deterministic_params,
      rng_stream = RNGStreams$new(),
      simulate_at_initialization = TRUE
    ) {
      #browser()
      private$.rng_stream <- rng_stream
      private$.rng_stream$resetAllStreams()
      private$.rng_stream$activateStream("population")
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )

      private$.target_population_size <- target_population_size
      private$.parameters <- parameters
      private$.generatePopulationFromRecipe(
        recipe = recipe,
        drop_extra_variables = FALSE,
        drop_empty_strata = TRUE,
        check = TRUE
      )
      private$.updatePersonParameterValues(
        reset_rng_stream = TRUE,
        drop_intermediate_params = TRUE
      )
      private$.simulation_erased <- TRUE

      if (simulate_at_initialization) {
        self$update()
      }
      invisible(self)
    },

    #' @description Run the natural-history simulation, filling event columns.
    update = function() {
      stopifnot(private$.simulation_erased)
      private$.generateAgeDeadFromOtherCauses(reset_rng_stream = TRUE)
      private$.generateDiseaseOnset(reset_rng_stream = TRUE)
      private$.generateDiseaseProgression(reset_rng_stream = TRUE)
      private$.generateAgeDeadFromLungCancer(reset_rng_stream = TRUE)
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

      private$.recipe[, `:=`(
        sim_count = ceiling(
          weight / sum(weight) * private$.target_population_size
        ),
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

    .setPopulation = function(population, copy) {
      stopifnot(data.table::is.data.table(population))
      required <- c("id", "sex", "spawn_age", "max_simulation_age")
      stopifnot(all(required %in% names(population)))
      private$.population <- if (copy) {
        data.table::copy(population)
      } else {
        population
      }
      data.table::setindex(private$.population, id)
      invisible(self)
    },

    .updatePersonParameterValues = function(
      reset_rng_stream,
      drop_intermediate_params = TRUE
    ) {
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )
      #browser()
      if (reset_rng_stream) {
        private$.rng_stream$resetStream("parameters")
      }
      private$.rng_stream$activateStream("parameters")

      private$.population <- private$.parameters[
        private$.population,
        on = c("sex", "smoking_status_at_spawn_age")
      ]

      private$.population[,
        multiplier := (1 + alpha1 * smoking_intensity_at_spawn_age^alpha2)
      ]
      private$.population[, cell_division_rate := alpha0 * multiplier]
      private$.population[, gamma := gamma0 * multiplier]
      private$.population[,
        cell_apoptosis_rate := cell_division_rate -
          malignant_cell_conversion_rate -
          gamma
      ]

      if (drop_intermediate_params) {
        drop_cols <- c(
          "multiplier",
          "alpha0",
          "alpha1",
          "alpha2",
          "gamma",
          "gamma0"
        )
        private$.population[, (drop_cols) := NULL]
      }

      invisible(self)
    },

    .generateAgeDeadFromOtherCauses = function(reset_rng_stream) {
      on.exit(
        private$.rng_stream$activateStream("population"),
        add = TRUE,
        after = FALSE
      )
      if (reset_rng_stream) {
        private$.rng_stream$resetStream("death_from_other_causes")
      }
      private$.rng_stream$activateStream("death_from_other_causes")

      #browser()

      private$.population[,
        age_dead_from_other_causes := generate_death_other(
          age_entry = spawn_age,
          sex = sex,
          smoking_status = smoking_status_at_spawn_age,
          intensity = smoking_intensity_at_spawn_age,
          tail = getOption("milc.death_other_tail")
        )
      ]

      if ("age_dead_from_lung_cancer" %in% names(private$.population)) {
        private$.population[,
          age_dead := pmin(
            age_dead_from_other_causes,
            age_dead_from_lung_cancer,
            na.rm = TRUE
          )
        ]
      } else {
        private$.population[, age_dead := age_dead_from_other_causes]
      }

      private$.simulated_cols <- c(
        private$.simulated_cols,
        "age_dead_from_other_causes"
      )
      invisible(self)
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

      #browser()
      private$.population[,
        age_at_disease_onset := runif(
          .N,
          min = spawn_age,
          max = age_dead_from_other_causes
        )
      ]

      private$.resolveAges(index_var = "age_at_disease_onset")

      private$.simulated_cols <- c(
        private$.simulated_cols,
        "age_at_disease_onset"
      )
      invisible(self)
    },

    .getAgeFromTumorVolume = function(tumor_logvolume_var, age_var) {
      private$.population[,
        (age_var) := age_at_disease_onset +
          (-tumor_growth_scale^(-1) *
            log(
              1 -
                (tumor_growth_scale / tumor_growth_shape) *
                  (get(tumor_logvolume_var) - tumor_logvolume_min)
            ))
      ]
      invisible(self)
    },

    .applyLogvolumeBounds = function(tumor_logvolume_var) {
      #browser()
      private$.population[,
        (tumor_logvolume_var) := fifelse(
          (get(tumor_logvolume_var) < tumor_logvolume_min) |
            (get(tumor_logvolume_var) > tumor_logvolume_max),
          NA,
          get(tumor_logvolume_var)
        )
      ]
      invisible(self)
    },

    # .applyAgeBoundsForTumorVolumes = function(age_var) {
    #   #browser()
    #   private$.population[,
    #     (age_var) := fifelse(
    #       (get(age_var) < spawn_age) |
    #         (get(age_var) > age_dead_from_other_causes),
    #       NA,
    #       get(age_var)
    #     )
    #   ]
    #   if (age_var != "age_at_disease_onset") {
    #     private$.population[,
    #       (age_var) := fifelse(
    #         (get(age_var) < age_at_disease_onset),
    #         NA,
    #         get(age_var)
    #       )
    #     ]
    #   }
    #   invisible(self)
    # },

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

      #browser()

      private$.population[,
        u_regional := runif(.N, min = 0, max = 1)
      ]
      private$.population[,
        u_distal := runif(.N, min = 0, max = 1)
      ]
      # sample logvolume_regional
      private$.population[,
        tumor_logvolume_regional := qnorm(
          u_regional,
          mean = tumor_logvolume_regional_mean,
          sd = tumor_logvolume_regional_sd
        )
      ]
      # sample logvolume_distal as conditional normal distribution
      # so that tumor_logvolume_distal > tumor_logvolume_regional
      private$.population[,
        tumor_logvolume_distal := qnorm(
          u_regional + u_distal * (1 - u_regional),
          mean = tumor_logvolume_distal_mean,
          sd = tumor_logvolume_distal_sd
        )
      ]
      # sample logvolume_diagnosis
      private$.population[,
        tumor_logvolume_diagnosis := rnorm(
          .N,
          mean = tumor_logvolume_diagnosis_mean,
          sd = tumor_logvolume_diagnosis_sd
        )
      ]

      # cleanup
      private$.applyLogvolumeBounds(
        tumor_logvolume_var = "tumor_logvolume_regional"
      )
      private$.applyLogvolumeBounds(
        tumor_logvolume_var = "tumor_logvolume_distal"
      )
      private$.applyLogvolumeBounds(
        tumor_logvolume_var = "tumor_logvolume_diagnosis"
      )

      private$.getAgeFromTumorVolume(
        tumor_logvolume_var = "tumor_logvolume_regional",
        age_var = "age_at_regional_disease"
      )
      private$.getAgeFromTumorVolume(
        tumor_logvolume_var = "tumor_logvolume_distal",
        age_var = "age_at_distal_disease"
      )
      private$.getAgeFromTumorVolume(
        tumor_logvolume_var = "tumor_logvolume_diagnosis",
        age_var = "age_at_diagnosis"
      )
      # private$.applyAgeBoundsForTumorVolumes(
      #   age_var = "age_at_regional_disease"
      # )
      # private$.applyAgeBoundsForTumorVolumes(age_var = "age_at_distal_disease")
      # private$.applyAgeBoundsForTumorVolumes(age_var = "age_at_diagnosis")
      private$.resolveAges(index_var = "age_at_regional_disease")
      private$.resolveAges(index_var = "age_at_distal_disease")
      private$.resolveAges(index_var = "age_at_diagnosis")

      private$.population[,
        stage_at_diagnosis := fcase(
          age_at_diagnosis < pmin(age_at_regional_disease, Inf, na.rm = TRUE)                                                                     , "local_disease"    ,
          age_at_diagnosis > pmin(age_at_regional_disease, Inf, na.rm = TRUE) & age_at_diagnosis < pmin(age_at_distal_disease, Inf, na.rm = TRUE) , "regional_disease" ,
          age_at_diagnosis > pmin(age_at_distal_disease, Inf, na.rm = TRUE)                                                                       , "distal_disease"   ,
          default = NA_character_
        )
      ]

      drop_cols <- c(
        "u_regional",
        "u_distal",
        "tumor_logvolume_regional",
        "tumor_logvolume_distal",
        "tumor_logvolume_diagnosis"
      )

      private$.population[, (drop_cols) := NULL]

      private$.simulated_cols <- c(
        private$.simulated_cols,
        "age_at_regional_disease",
        "age_at_distal_disease",
        "age_at_diagnosis",
        "stage_at_diagnosis"
      )
      invisible(self)
    },

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

      private$.population[, age_dead_from_lung_cancer := NA_real_]
      private$.population[
        !is.na(age_at_disease_onset),
        age_dead_from_lung_cancer := runif(
          n = .N,
          min = age_at_disease_onset,
          max = getOption("milc.max_age")
        )
      ]
      private$.population[,
        age_dead := pmin(
          age_dead_from_other_causes,
          age_dead_from_lung_cancer,
          na.rm = TRUE
        )
      ]

      ### now resolve all ages
      private$.resolveAges(index_var = "age_at_disease_onset")
      private$.resolveAges(index_var = "age_at_diagnosis")
      private$.resolveAges(index_var = "age_at_regional_disease")
      private$.resolveAges(index_var = "age_at_distal_disease")
      private$.resolveAges(index_var = "age_dead")
      private$.resolveAges(index_var = "age_dead_from_lung_cancer")

      private$.population[,
        cause_of_death := fifelse(
          is.na(age_dead_from_lung_cancer),
          "other_causes",
          "lung_cancer"
        )
      ]

      private$.simulated_cols <- c(
        private$.simulated_cols,
        "age_dead_from_lung_cancer",
        "cause_of_death"
      )
      invisible(self)
    },

    .resolveAges = function(index_var) {
      stopifnot(!is.null(index_var))
      disease_progression_vars <- c(
        "age_at_regional_disease",
        "age_at_distal_disease",
        "age_at_diagnosis",
        "age_dead_from_lung_cancer",
        "age_dead"
      )
      all_age_vars <- c(
        disease_progression_vars,
        "age_at_disease_onset",
        "age_dead_from_other_causes"
      )
      existing_vars <- intersect(
        x = all_age_vars,
        y = names(private$.population)
      )
      stopifnot(index_var %in% existing_vars)
      missing_vars <- setdiff(x = all_age_vars, y = names(private$.population))

      private$.population[,
        (index_var) := fifelse(
          get(index_var) %between% .(spawn_age, getOption("milc.max_age")),
          get(index_var),
          NA
        )
      ]

      if (index_var %in% disease_progression_vars) {
        stopifnot("age_at_disease_onset" %in% existing_vars)
        stopifnot("age_dead" %in% existing_vars)
        switch(
          EXPR = index_var,
          age_at_diagnosis = ,
          age_at_regional_disease = {
            private$.population[,
              (index_var) := fifelse(
                get(index_var) %between%
                  .(age_at_disease_onset, age_dead),
                get(index_var),
                NA
              )
            ]
          },
          age_at_distal_disease = {
            stopifnot("age_at_regional_disease" %in% existing_vars)
            private$.population[,
              (index_var) := fifelse(
                get(index_var) %between%
                  .(age_at_regional_disease, age_dead),
                get(index_var),
                NA
              )
            ]
          },
          age_dead_from_lung_cancer = {
            stopifnot("age_dead_from_other_causes" %in% existing_vars)
            private$.population[,
              (index_var) := fifelse(
                get(index_var) %between%
                  .(age_at_regional_disease, age_dead_from_other_causes),
                get(index_var),
                NA
              )
            ]
          }
        )
      }

      invisible(self)
    }
  )
)
