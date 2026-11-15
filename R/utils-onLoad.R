#' @import data.table
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.milc <- list(
    milc.verbose = FALSE,
    milc.eps = .Machine$double.eps,
    milc.exp_max = log(.Machine$double.xmax) - log(10^15),

    milc.max_age = 100, # years
    milc.death_other_tail = "archival", # or "corrected"; see generate_death_other()
    milc.progression_method = "rejection", # or "conditional" (CRN); see generate_progression()

    milc.globalRNG.seed = NULL,
    milc.antithetic = FALSE,
    milc.rng_stream_names = c(
      "simulation",
      "population",
      "parameters",
      "death_from_other_causes",
      "disease_onset",
      "disease_progression",
      "death_from_lung_cancer",
      "screening",
      "treatment",
      "calibration_design_generation"
    )
  )

  toset <- !(names(op.milc) %in% names(op))
  if (any(toset)) {
    options(op.milc[toset])
  }
  remove_lecuyer()
  invisible()
}


.onUnload <- function(libname, pkgname) {
  remove_lecuyer()
  invisible()
}
