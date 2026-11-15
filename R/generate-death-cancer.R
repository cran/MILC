#' Map stage at diagnosis to the archival labels
#'
#' Accepts the archival labels `"localized"`, `"regional"`, `"distant"` and the
#' labels used by the `Population` class, `"local_disease"`,
#' `"regional_disease"`, `"distal_disease"`.
#'
#' @param stage character vector.
#' @return character vector in `c("localized", "regional", "distant")`.
#' @keywords internal
milc_stage <- function(stage) {
  map <- c(
    localized = "localized", regional = "regional", distant = "distant",
    local_disease = "localized", regional_disease = "regional",
    distal_disease = "distant"
  )
  out <- unname(map[stage])
  bad <- unique(stage[is.na(out)])
  if (length(bad) > 0L) {
    stop("Unknown stage: ", paste(bad, collapse = ", "))
  }
  out
}

#' Row of `rates_death_lung` for each individual
#'
#' The stratum is (stage at diagnosis, age group at diagnosis, sex). Where the
#' archival `age_grp` returns `NA` (age exactly 120), the literal `"NA"` fallback
#' strata are used, as the archival `grep` lookup does; only `"distant NA male"`
#' and `"localized NA female"` exist.
#'
#' @inheritParams generate_death_cancer
#' @return integer vector of rows of `rates$lambda_matrix`.
#' @keywords internal
death_cancer_stratum <- function(age_diagnosis, stage, sex,
                                 rates = fastMILC::rates_death_lung) {
  agegrp <- milc_age_group(age_diagnosis)
  agegrp[is.na(agegrp)] <- "NA"
  keys <- data.table(
    stage = milc_stage(stage),
    agegrp = agegrp,
    sex = tolower(sex)
  )
  row <- rates$strata[keys, on = c("stage", "agegrp", "sex"), x.stratum_row]
  if (anyNA(row)) {
    miss <- unique(keys[is.na(row)])
    stop(
      "No lung-cancer mortality table for: ",
      paste(do.call(paste, miss), collapse = "; ")
    )
  }
  row
}

#' Generate age at death from lung cancer
#'
#' Samples the latent age at death from lung cancer after clinical diagnosis, as
#' in the archival `tdeath_lung` (Chrysanthopoulou, 2014), using the step rates
#' in [rates_death_lung] and [nhppp::vdraw_sc_step()].
#'
#' Lung-cancer death is possible after diagnosis at any stage (archival; IJM
#' instead requires the distant state). The clock is time since diagnosis and
#' the stratum is fixed by stage, age group at diagnosis, and sex. Within the
#' tabulated window the first arrival reproduces the archival cumulative
#' incidence at every positive knot exactly; the atom at diagnosis is merged
#' into the first cell. If there is no event in the window:
#' \itemize{
#'   \item in a stratum whose cumulative incidence reaches 1
#'     (`terminal_certain`), the time is uniform in `[terminal_lo, horizon]`;
#'   \item otherwise (probability `1 - F_max`), death from lung cancer occurs at
#'     `age_diagnosis + horizon`, the archival point mass (`tdeath_lung` returns
#'     the horizon with `event = 0`, which `nat_hist` does not use).
#' }
#' The stratum `"localized NA female"` has no cells (horizon 0), so its draws
#' equal `age_diagnosis`, as in the archival code.
#'
#' @param age_diagnosis numeric vector, age at clinical diagnosis (years); `NA`
#'   for individuals never diagnosed.
#' @param stage character vector, stage at diagnosis; see [milc_stage()] for
#'   accepted labels. Ignored where `age_diagnosis` is `NA`.
#' @param sex character vector, `"male"` or `"female"`.
#' @param rates the rate table, default [rates_death_lung].
#' @param chunk_size maximum number of rows per [nhppp::vdraw_sc_step()] call.
#' @return numeric vector, age at death from lung cancer; `NA` where
#'   `age_diagnosis` is `NA`.
#' @export
generate_death_cancer <- function(age_diagnosis, stage, sex,
                                  rates = fastMILC::rates_death_lung,
                                  chunk_size = 1e5L) {
  n <- length(age_diagnosis)
  stopifnot(
    is.numeric(age_diagnosis), length(stage) == n, length(sex) == n,
    chunk_size >= 1
  )
  chunk_size <- as.integer(chunk_size)
  out <- rep(NA_real_, n)
  dx <- which(!is.na(age_diagnosis))
  if (length(dx) == 0L) {
    return(out)
  }
  if (anyNA(stage[dx])) {
    stop("`stage` is NA for ", sum(is.na(stage[dx])), " diagnosed individual(s)")
  }
  row <- death_cancer_stratum(age_diagnosis[dx], stage[dx], sex[dx], rates)

  # Strata in increasing row order: deterministic random-number consumption.
  tau <- rep(NA_real_, length(dx))
  by_row <- split(seq_along(dx), row)
  for (r in names(by_row)) {
    idx <- by_row[[r]]
    tau[idx] <- draw_step_first_arrival(
      rates, as.integer(r), length(idx), chunk_size
    )
  }

  st <- rates$strata[row]
  none <- is.na(tau)
  certain <- none & st$terminal_certain
  if (any(certain)) {
    tau[certain] <- stats::runif(
      sum(certain), st$terminal_lo[certain], st$horizon[certain]
    )
  }
  residual <- none & !st$terminal_certain
  tau[residual] <- st$horizon[residual]

  out[dx] <- age_diagnosis[dx] + tau
  out
}
