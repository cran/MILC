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

# Progression (archival t_prog parameterization; model.tex, Section 5).
# Values: IJM Table 2 posterior medians (m, mu_reg, mu_dist, mu_diagn), with
# location = scale (IJM). The archival code has no defaults (the nat_hist
# examples use m = 0.00042, cdiagn = 3.91, creg = 1.1, cdist = 2.8).
#   tumor_growth_rate          m        Gompertz scale, per day
#   volume_diagnosis_meanlog   cdiagn[1], volume_diagnosis_sdlog cdiagn[2]
#   volume_regional_meanlog    creg[1],   volume_regional_sdlog  creg[2]
#   volume_distant_meanlog     cdist[1],  volume_distant_sdlog   cdist[2]
# (log-normal parameters of volume in ml). generate_progression() reproduces
# the archival swap: V_reg is drawn with cdist and V_dist with creg.
deterministic_params[, `:=`(
  tumor_growth_rate = 3.4e-4,
  volume_diagnosis_meanlog = 2.65,
  volume_diagnosis_sdlog = 2.65,
  volume_regional_meanlog = 2.16,
  volume_regional_sdlog = 2.16,
  volume_distant_meanlog = 5.62,
  volume_distant_sdlog = 5.62
)]

usethis::use_data(deterministic_params, overwrite = TRUE)
