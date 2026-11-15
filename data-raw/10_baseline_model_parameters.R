# Baseline model parameters -> package dataset `deterministic_params`.
#
# One row per sex. Onset (TSCE) parameters are the archival values of
# archival_code/R/HT_mal.R, which fastMILC follows where they differ from
# Chrysanthopoulou, Int J Microsimulation 2017, Table 2 (model.tex, Table 1
# and Section 9). Symbols as in model.tex, Section 4 (archival names in
# brackets):
#   n_stem_cells                    X        [X = 10^7]
#   normal_cell_initiation_rate     nu_ns    [v0]
#   malignant_cell_conversion_rate  mu       [m = v0]
#   cell_division_rate_nonsmoking   alpha_ns [a0, smok_status "N"]
#   cell_division_rate_smoking_base alpha_s0 [a0, smok_status "Y"]
#   net_proliferation_rate          gamma    [g0; gamma = alpha - beta - mu]
#   smoking_initiation_multiplier   rho      [1 + p4]
#   smoking_dose_coefficient        alpha_1  [p5]
#   smoking_dose_exponent           alpha_2  [p6]
# Smoking parameter set: kappa(q) = 1 + alpha_1 q^alpha_2, alpha_s =
# alpha_s0 kappa(q), gamma_s = gamma kappa(q), nu_s = rho nu_ns, mu_s = mu.
#
# Progression parameters are unchanged from the previous version (IJM values,
# log-volumes for diameters in cm); step 5 replaces them with the archival
# t_prog parameterization.
#
# Run from the package root: Rscript data-raw/10_baseline_model_parameters.R

library(data.table)

deterministic_params <- data.table(
  sex = c("male", "female"),
  n_stem_cells = 1e7,
  normal_cell_initiation_rate = c(7e-8, 1e-7),
  malignant_cell_conversion_rate = c(7e-8, 1e-7),
  cell_division_rate_nonsmoking = c(7.7, 15.5),
  cell_division_rate_smoking_base = c(7.5, 15.5),
  net_proliferation_rate = c(0.09, 0.071),
  smoking_initiation_multiplier = c(1.01, 1.02),
  smoking_dose_coefficient = c(0.6, 0.5),
  smoking_dose_exponent = c(0.22, 0.32)
)

deterministic_params[, `:=`(
  tumor_logvolume_min = log(pi / 6 * (10^-5 * 100)^3), # diameters in cm!
  tumor_logvolume_max = log(pi / 6 * (0.13 * 100)^3), # diameters in cm!

  tumor_growth_scale = (3.4 * 10^-4) * 365.24, # scale time_intervals in YEARS!
  tumor_growth_shape = 31 * (3.4 * 10^-4) * 365.24, # scale time_intervals in YEARS!

  tumor_logvolume_regional_mean = 2.16, # for diameters in cm!
  tumor_logvolume_regional_sd = 2.16, # for diameters in cm!

  tumor_logvolume_distal_mean = 5.62, # for diameters in cm!
  tumor_logvolume_distal_sd = 5.62, # for diameters in cm!

  tumor_logvolume_diagnosis_mean = 2.65, # for diameters in cm!
  tumor_logvolume_diagnosis_sd = 2.65 # for diameters in cm!
)]

usethis::use_data(deterministic_params, overwrite = TRUE)
