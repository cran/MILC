# Tumor growth and disease progression (model_description/model.tex,
# Section 5): Gompertz growth from one malignant cell and three random volume
# thresholds (regional, distant, diagnosis), drawn exactly from their joint
# truncated log-normal law.

# Volume of one malignant cell and the maximum volume (mm^3), archival t_prog.
MILC_V0 <- pi * 0.01^3 / 6
MILC_VMAX <- pi * 130^3 / 6

#' Time from onset to a tumor volume (Gompertz growth)
#'
#' Vectorized archival `t_vol`: \eqn{\tau(v) = -(365 m)^{-1}\log(1 -
#' \log(v/V_0)/31)} years, with \eqn{V_0 = \pi 0.01^3/6} mm\eqn{^3} and shape
#' \eqn{s = 31m}.
#'
#' @param vol volume (mm\eqn{^3}).
#' @param m Gompertz scale parameter (per day).
#' @return years from onset.
#' @keywords internal
milc_time_to_volume <- function(vol, m) {
  -log1p(-log(vol / MILC_V0) / 31) / (365 * m)
}

#' Draw volumes from a log-normal truncated to \eqn{(V_0, V_{max})}
#'
#' Volume \eqn{V = 1000 Z} mm\eqn{^3}, \eqn{Z \sim LN(\mu, \sigma)} in ml,
#' truncated to \eqn{V_0 < V < V_{max}}, by inverse CDF (exact).
#'
#' @param meanlog,sdlog log-normal parameters of \eqn{Z}, one per draw.
#' @return volumes (mm\eqn{^3}).
#' @keywords internal
draw_truncated_volume <- function(meanlog, sdlog) {
  n <- length(meanlog)
  lo <- stats::pnorm((log(MILC_V0 / 1000) - meanlog) / sdlog)
  hi <- stats::pnorm((log(MILC_VMAX / 1000) - meanlog) / sdlog)
  p <- lo + stats::runif(n) * (hi - lo)
  v <- 1000 * exp(meanlog + sdlog * stats::qnorm(p))
  # guard against rounding onto the bounds
  pmin(pmax(v, MILC_V0 * (1 + 1e-12)), MILC_VMAX * (1 - 1e-12))
}

#' Generate disease progression after onset
#'
#' Samples the volume thresholds \eqn{(V_{reg}, V_{dist}, V_{diagn})} from the
#' joint law of independent log-normals truncated to
#' \eqn{\{V_0 < V_{diagn} < V_{max}\} \cap \{V_0 < V_{reg} < V_{dist} <
#' V_{max}\}} (`model.tex`, Eq. truncset): \eqn{V_{diagn}} by inverse CDF of
#' its truncated law; \eqn{V_{reg}} and \eqn{V_{dist}} by inverse CDF of their
#' box-truncated laws, redrawing only the rows with \eqn{V_{reg} \ge
#' V_{dist}}. This is exact, as is the archival whole-vector rejection in
#' `t_prog`.
#'
#' As in the archival `t_prog`, \eqn{V_{reg}} is drawn with the parameters
#' labelled `distant` (`cdist`) and \eqn{V_{dist}} with those labelled
#' `regional` (`creg`); see `model.tex`, Sections 5 and 9.
#'
#' @param age_onset age at onset of the first malignant cell; `NA` for no onset.
#' @param params data.frame with one row per person and columns
#'   `tumor_growth_rate` (\eqn{m}, per day), `volume_diagnosis_meanlog`,
#'   `volume_diagnosis_sdlog`, `volume_regional_meanlog`,
#'   `volume_regional_sdlog`, `volume_distant_meanlog`, `volume_distant_sdlog`
#'   (log-normal parameters of the volume in ml; archival `cdiagn`, `creg`,
#'   `cdist`).
#' @param max_rounds maximum rounds of redrawing for the order constraint.
#' @return data.table with one row per person: `age_at_regional_disease`,
#'   `age_at_distal_disease`, `age_at_diagnosis`, `diameter_at_diagnosis` (mm),
#'   `stage_at_diagnosis` (`"localized"`, `"regional"`, `"distant"`); all `NA`
#'   where `age_onset` is `NA`. These are latent ages: diagnosis may fall after
#'   death.
#' @export
generate_progression <- function(age_onset, params, max_rounds = 1000L) {
  n <- length(age_onset)
  stopifnot(nrow(params) == n)
  out <- data.table(
    age_at_regional_disease = rep(NA_real_, n),
    age_at_distal_disease = rep(NA_real_, n),
    age_at_diagnosis = rep(NA_real_, n),
    diameter_at_diagnosis = rep(NA_real_, n),
    stage_at_diagnosis = rep(NA_character_, n)
  )
  idx <- which(!is.na(age_onset))
  if (length(idx) == 0L) {
    return(out)
  }
  p <- params[idx, ]
  k <- length(idx)

  v_diagn <- draw_truncated_volume(
    p$volume_diagnosis_meanlog, p$volume_diagnosis_sdlog
  )
  # archival swap: V_reg uses cdist, V_dist uses creg
  v_reg <- rep(NA_real_, k)
  v_dist <- rep(NA_real_, k)
  todo <- seq_len(k)
  for (round in seq_len(max_rounds)) {
    v_reg[todo] <- draw_truncated_volume(
      p$volume_distant_meanlog[todo], p$volume_distant_sdlog[todo]
    )
    v_dist[todo] <- draw_truncated_volume(
      p$volume_regional_meanlog[todo], p$volume_regional_sdlog[todo]
    )
    todo <- todo[v_reg[todo] >= v_dist[todo]]
    if (length(todo) == 0L) break
  }
  if (length(todo) > 0L) {
    stop("Order constraint V_reg < V_dist not met after ", max_rounds,
         " rounds for ", length(todo), " individual(s)")
  }

  m <- p$tumor_growth_rate
  a <- age_onset[idx]
  set(out, idx, "age_at_regional_disease", a + milc_time_to_volume(v_reg, m))
  set(out, idx, "age_at_distal_disease", a + milc_time_to_volume(v_dist, m))
  set(out, idx, "age_at_diagnosis", a + milc_time_to_volume(v_diagn, m))
  set(out, idx, "diameter_at_diagnosis", (6 * v_diagn / pi)^(1 / 3))
  set(out, idx, "stage_at_diagnosis", fifelse(
    v_diagn < v_reg, "localized",
    fifelse(v_diagn < v_dist, "regional", "distant")
  ))
  out
}
