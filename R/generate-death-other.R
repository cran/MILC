#' Archival age-group labels (vectorized `age_grp`)
#'
#' Labels `"<40"`, `"40-44"`, ..., `"80-84"`, `">85"` with left-closed bands
#' (`ag(x)` in `model_description/model.tex`, Section 3). As in the archival
#' `age_grp`, age exactly 120 maps to `NA`, and ages above 120 to `">85"`.
#'
#' @param age numeric vector of ages (years).
#' @return character vector of archival age-group labels.
#' @keywords internal
milc_age_group <- function(age) {
  labels <- c(
    "<40", "40-44", "45-49", "50-54", "55-59", "60-64",
    "65-69", "70-74", "75-79", "80-84", ">85"
  )
  out <- labels[findInterval(age, seq(40, 85, by = 5)) + 1L]
  out[is.na(age) | age == 120] <- NA_character_
  out
}

#' Archival smoking-intensity labels (vectorized `d_grp`)
#'
#' `"NA"` for `q < 1`, `"1-10"` for `1 <= q < 11`, `"11-20"` for
#' `11 <= q < 21`, `"21+"` for `q >= 21` (cigarettes per day; `dg(q)` in
#' `model_description/model.tex`, Section 3).
#'
#' @param q numeric vector of cigarettes per day.
#' @return character vector of archival intensity labels (`NA` where `q` is `NA`).
#' @keywords internal
milc_intensity_group <- function(q) {
  labels <- c("NA", "1-10", "11-20", "21+")
  labels[findInterval(q, c(1, 11, 21)) + 1L]
}

#' Map smoking status to the archival labels
#'
#' Accepts the archival labels `"never"`, `"former"`, `"current"`, and the
#' recipe labels `"non_smoker"` (mapped to `"never"`) and `"smoker"` (mapped to
#' `"current"`).
#'
#' @param status character vector.
#' @return character vector in `c("never", "former", "current")`.
#' @keywords internal
milc_smoking_status <- function(status) {
  map <- c(
    never = "never", former = "former", current = "current",
    non_smoker = "never", smoker = "current"
  )
  out <- unname(map[status])
  bad <- unique(status[is.na(out)])
  if (length(bad) > 0L) {
    stop("Unknown smoking status: ", paste(bad, collapse = ", "))
  }
  out
}

#' Row of `rates_death_other` for each individual
#'
#' The stratum is (smoking status, age group at entry, sex, intensity group);
#' the intensity group is used for current smokers only (archival).
#'
#' @inheritParams generate_death_other
#' @return integer vector of rows of `rates$lambda_matrix`.
#' @keywords internal
death_other_stratum <- function(age_entry, sex, smoking_status, intensity,
                                rates = fastMILC::rates_death_other) {
  status <- milc_smoking_status(smoking_status)
  keys <- data.table(
    smoking_status = status,
    agegrp = milc_age_group(age_entry),
    sex = tolower(sex),
    intensity = fifelse(
      status == "current", milc_intensity_group(intensity), NA_character_
    )
  )
  row <- rates$strata[
    keys,
    on = c("smoking_status", "agegrp", "sex", "intensity"),
    x.stratum_row
  ]
  if (anyNA(row)) {
    miss <- unique(keys[is.na(row)])
    stop(
      "No other-cause mortality table for: ",
      paste(do.call(paste, miss), collapse = "; ")
    )
  }
  row
}

#' Sample latent first arrivals for individuals that share one rate row
#'
#' Calls [nhppp::vdraw_sc_step()] on the cells up to the row's horizon, in
#' chunks of at most `chunk_size` rows, so the rate matrix passed to nhppp is at
#' most `chunk_size` x (cells up to the horizon).
#'
#' @param rates a rate table (`rates_death_other` or `rates_death_lung`).
#' @param row row of `rates$lambda_matrix`.
#' @param n number of draws.
#' @param chunk_size maximum rows per nhppp call.
#' @return numeric vector of length `n`: first-arrival times on the table clock
#'   (years), `NA` when there is no event within the horizon.
#' @keywords internal
draw_step_first_arrival <- function(rates, row, n, chunk_size) {
  out <- rep(NA_real_, n)
  tb <- rates$time_breaks
  K <- which.min(abs(tb - rates$strata$horizon[row])) - 1L
  if (n == 0L || K == 0L) {
    return(out)
  }
  lam <- rates$lambda_matrix[row, seq_len(K)]
  if (all(lam == 0)) {
    return(out)
  }
  tb <- tb[seq_len(K + 1L)]
  starts <- seq.int(1L, n, by = chunk_size)
  for (s in starts) {
    idx <- s:min(n, s + chunk_size - 1L)
    out[idx] <- nhppp::vdraw_sc_step(
      lambda_matrix = matrix(lam, nrow = length(idx), ncol = K, byrow = TRUE),
      time_breaks = tb,
      atmost1 = TRUE
    )[, 1L]
  }
  out
}

#' Generate age at death from other causes
#'
#' Samples the latent age at death from causes other than lung cancer, as in
#' the archival `tdeath_other` (Chrysanthopoulou, 2014), using the step rates
#' in [rates_death_other] and [nhppp::vdraw_sc_step()].
#'
#' Within the stratum's tabulated window `[0, horizon]` (years since entry),
#' the first arrival of the Poisson process reproduces the archival cumulative
#' incidence at every knot exactly. Individuals in a stratum whose cumulative
#' incidence reaches 1 (`terminal_certain`) and who have no earlier event get a
#' uniform time in the terminal cell `[terminal_lo, horizon]`. Otherwise,
#' individuals with no event in the window (probability `1 - F_max`) get a
#' time from the tail rule:
#' \describe{
#'   \item{`"archival"` (default)}{`age_entry + U * (max_age - age_entry -
#'     horizon)`, i.e. `Unif(age_entry, max_age - horizon)`, the formula of the
#'     archival `tdeath_other`. This overlaps the tabulated window. For
#'     `age_entry > max_age - horizon` the width is negative and the age at
#'     death precedes entry, as in the archival code (whose `nat_hist` then
#'     excludes the record).}
#'   \item{`"corrected"`}{`Unif(age_entry + horizon, max_age)`, the presumably
#'     intended rule. For `age_entry + horizon > max_age` the draw exceeds
#'     `max_age` (with a warning).}
#' }
#' The archival code uses `max_age = 110`; fastMILC uses the package option
#' `milc.max_age` (100). The model parameters were calibrated under the archival
#' rule (`model_description/model.tex`, Sections 2 and 6).
#'
#' @param age_entry numeric vector, age at entry (years). Fixes the age group of
#'   the table for the whole follow-up (archival).
#' @param sex character vector, `"male"` or `"female"`.
#' @param smoking_status character vector, smoking status at entry; see
#'   [milc_smoking_status()] for accepted labels.
#' @param intensity numeric vector, cigarettes per day at entry; used for
#'   current smokers only.
#' @param tail `"archival"` or `"corrected"`.
#' @param max_age upper age of the tail rule; default `getOption("milc.max_age")`
#'   (100). The archival code uses 110.
#' @param rates the rate table, default [rates_death_other].
#' @param chunk_size maximum number of rows per [nhppp::vdraw_sc_step()] call.
#' @return numeric vector, age at death from other causes.
#' @export
generate_death_other <- function(age_entry, sex, smoking_status, intensity,
                                 tail = c("archival", "corrected"),
                                 max_age = getOption("milc.max_age"),
                                 rates = fastMILC::rates_death_other,
                                 chunk_size = 1e5L) {
  tail <- match.arg(tail)
  n <- length(age_entry)
  stopifnot(
    is.numeric(age_entry), !anyNA(age_entry),
    length(sex) == n, length(smoking_status) == n, length(intensity) == n,
    is.numeric(max_age), length(max_age) == 1L,
    chunk_size >= 1
  )
  chunk_size <- as.integer(chunk_size)
  row <- death_other_stratum(age_entry, sex, smoking_status, intensity, rates)

  # Duration since entry. Strata are processed in increasing row order, so the
  # random-number consumption is deterministic given the population.
  tau <- rep(NA_real_, n)
  by_row <- split(seq_len(n), row)
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

  age <- age_entry + tau
  in_tail <- none & !st$terminal_certain
  if (any(in_tail)) {
    a <- age_entry[in_tail]
    h <- st$horizon[in_tail]
    u <- stats::runif(sum(in_tail))
    if (tail == "corrected" && any(a + h > max_age)) {
      warning(
        sum(a + h > max_age), " tail draw(s) with age_entry + horizon > max_age; ",
        "the corrected tail Unif(age_entry + horizon, max_age) is then reversed"
      )
    }
    age[in_tail] <- switch(
      tail,
      archival = a + u * (max_age - (a + h)),
      corrected = (a + h) + u * (max_age - (a + h))
    )
  }
  age
}
