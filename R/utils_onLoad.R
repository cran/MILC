
.onLoad <- function(libname, pkgname){
  op <- options()
  op.pneumon <- list(
    pneumon.verbose = FALSE,
    pneumon.eps = .Machine$double.eps,
    pneumon.exp_max = log(.Machine$double.xmax) - log(10^15),
    pneumon.n_stem_cells = 10^7,          # cells  
    pneumon.max_age = 100,                # years
    pneumon.cell_diameter = 0.01,         # mm
    pneumon.maximum_tumor_diameter = 130, # mm     
    pneumon.globalRNG.seed = NULL,
    pneumon.antithetic = FALSE, 
    pneumon.rng_stream_names = c(
      "simulation", 
      "population", 
      "disease_onset", 
      "disease_evolution"
    )
  )


  toset <- !(names(op.pneumon) %in% names(op))
  if (any(toset)) options(op.pneumon[toset])
  remove_lecuyer()
  invisible()
}




.onUnload <- function(libname, pkgname) {
  remove_lecuyer()
  invisible()
}
