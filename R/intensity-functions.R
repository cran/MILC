# Ported, vectorized MILC kernel functions.
# These are faithful vectorizations of the archival scalar functions
# (archival_code/R/{HT_mal,ht_mal_int,t_vol,age_grp,d_grp}.R). They are pure and
# side-effect free; correctness is checked against the archival scalars in
# tests/testthat/test-intensity-functions.R.

# Two-stage clonal-expansion parameters by gender x smoking status.
# Rows: a0, g0, v0, p4, p5, p6  (columns keyed "<gender> <status>", status N/Y).
.milc_tsce_params <- local({
  m <- cbind(
    "male N"   = c(a0 = 7.7,  g0 = 0.090, v0 = 7e-8, p4 = 0.00, p5 = 0.0, p6 = 0.00),
    "female N" = c(a0 = 15.5, g0 = 0.071, v0 = 1e-7, p4 = 0.00, p5 = 0.0, p6 = 0.00),
    "male Y"   = c(a0 = 7.5,  g0 = 0.090, v0 = 7e-8, p4 = 0.01, p5 = 0.6, p6 = 0.22),
    "female Y" = c(a0 = 15.5, g0 = 0.071, v0 = 1e-7, p4 = 0.02, p5 = 0.5, p6 = 0.32)
  )
  m
})

#' MILC malignancy-onset cumulative intensity
#'
#' Closed-form cumulative intensity (integrated hazard) of the two-stage
#' clonal-expansion model for malignant-cell onset, per Chrysanthopoulou's MILC.
#' This is the \eqn{\Lambda(t)} that the onset sampler inverts via
#' [nhppp::vdraw_cumulative_intensity()]. Vectorized over all arguments (scalars
#' are recycled).
#'
#' @param t Age (years); scalar or vector.
#' @param gender `"male"` or `"female"`; scalar or vector.
#' @param d Smoking intensity (cigarettes/day). Unused for non-smokers; may be `NA`.
#' @param smok_status `"N"` (not currently accruing smoking risk) or `"Y"`.
#' @return Cumulative intensity \eqn{\Lambda(t)} (same length as the recycled inputs).
#' @export
milc_HT_mal <- function(t, gender, d, smok_status) {
  X <- 1e7
  key <- paste(gender, smok_status)
  if (!all(key %in% colnames(.milc_tsce_params))) {
    stop("Unknown gender/smok_status combination(s): ",
         paste(unique(key[!key %in% colnames(.milc_tsce_params)]), collapse = ", "))
  }
  P <- .milc_tsce_params[, key, drop = FALSE]     # 6 x n
  a0 <- P["a0", ]; g0 <- P["g0", ]; v0 <- P["v0", ]
  p4 <- P["p4", ]; p5 <- P["p5", ]; p6 <- P["p6", ]

  mult <- 1 + p5 * d^p6                           # NA^0 == 1, and p5==0 for non-smokers
  v <- v0 * (1 + p4)
  g <- g0 * mult
  a <- a0 * mult
  m <- v0
  B <- 0.5 * (-g + sqrt(g^2 + 4 * a * m))
  K <- v * m * X / (g + B)
  unname(K * (-t + (1 / B) * log(g + B + B * exp((2 * B + g) * t))))
}

#' MILC onset cumulative intensity over an interval
#'
#' Integrated onset intensity over `[lower, upper]`, i.e.
#' `milc_HT_mal(upper, ...) - milc_HT_mal(lower, ...)`. Vectorized.
#'
#' @inheritParams milc_HT_mal
#' @param lower,upper Interval endpoints (ages, years).
#' @return Integrated intensity over the interval.
#' @export
milc_ht_mal_int <- function(lower, upper, gender, d, smok_status) {
  milc_HT_mal(upper, gender, d, smok_status) -
    milc_HT_mal(lower, gender, d, smok_status)
}

#' Time for a tumor to reach a given volume (Gompertz growth inverse)
#'
#' Inverts the Gompertz tumor-growth model to return the time (years) for a
#' tumor to grow from the seed volume to `vol`. Vectorized over `vol` and `m`.
#'
#' @param vol Tumor volume (cubic mm).
#' @param m Gompertz scale parameter.
#' @return Time in years.
#' @export
milc_t_vol <- function(vol, m) {
  d0 <- 0.01
  V0 <- pi * (d0^3) / 6
  s  <- 31 * m
  (-1 / m) * log(1 - (m / s) * log(vol / V0)) / 365
}

#' MILC age group
#'
#' Maps age (years) to the MILC 5-year age band. Vectorized.
#'
#' @param age Age in years; scalar or vector.
#' @return Character age-group label (`"<40"`, `"40-44"`, ..., `">85"`), `NA` for `NA` age.
#' @export
milc_age_grp <- function(age) {
  age_char <- c("<40", "40-44", "45-49", "50-54", "55-59", "60-64",
                "65-69", "70-74", "75-79", "80-84", ">85")
  age_num  <- c(40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 120)
  idx <- findInterval(age, age_num) + 1L
  idx[idx > length(age_char)] <- length(age_char)   # age >= 120 -> ">85"
  out <- age_char[idx]
  out[is.na(age)] <- NA_character_
  out
}

#' MILC smoking-intensity group
#'
#' Maps smoking intensity (cigarettes/day) to the MILC intensity band. Vectorized.
#'
#' @param d Smoking intensity; scalar or vector.
#' @return Character label (`"1-10"`, `"11-20"`, `"21+"`; `"NA"` band for `d < 1`),
#'   `NA` for `NA` input.
#' @export
milc_d_grp <- function(d) {
  d_char <- c("NA", "1-10", "11-20", "21+")
  d_num  <- c(1, 11, 21)
  idx <- findInterval(d, d_num) + 1L
  out <- d_char[idx]
  out[is.na(d)] <- NA_character_
  out
}
