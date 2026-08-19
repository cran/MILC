# baseline model parameters
library(data.table)

# From Chrysanthopoulou Int J Microsimulation 2017, Table 2
deterministic_tumor_onset_male_nonsmoker <- data.table(
  sex = "male",
  smoking_status_at_spawn_age = "non_smoker",
  normal_cell_initiation_rate = 7.16 * 10^-8, # normal initiation rate
  alpha0 = 7.7, # division rate
  alpha1 = 0,
  alpha2 = 0,
  gamma0 = 0.09,
  malignant_cell_conversion_rate = 7.16 * 10^-8 # malignant_conversion_rate
)

deterministic_tumor_onset_female_nonsmoker <- data.table(
  sex = "female",
  smoking_status = "non_smoker",
  normal_cell_initiation_rate = 1.07 * 10^-7, # normal initiation rate
  alpha0 = 15.82, # division rate
  alpha1 = 0,
  alpha2 = 0,
  gamma0 = 0.071,
  malignant_cell_conversion_rate = 1.07 * 10^-7 # malignant_conversion_rate
)

deterministic_tumor_onset_male_smoker <- deterministic_params_male_nonsmoker <- data.table(
  sex = "male",
  smoking_status = "smoker",
  normal_cell_initiation_rate = 7.16 * 10^-8, # nu
  alpha0 = 7.7,
  alpha1 = 0.6,
  alpha2 = 0.22,
  gamma0 = 0.09,
  malignant_cell_conversion_rate = 7.16 * 10^-8 # mu
)

deterministic_tumor_onset_female_smoker <- data.table(
  sex = "female",
  smoking_status = "smoker",
  normal_cell_initiation_rate = 0.98 * 1.07 * 10^-7, # nu
  alpha0 = 15.82, # division rate
  alpha1 = 0.5,
  alpha2 = 0.32,
  gamma0 = 0.071,
  malignant_cell_conversion_rate = 1.07 * 10^-7 # mu
)

deterministic_params <- rbindlist(
  list(
    deterministic_tumor_onset_male_nonsmoker,
    deterministic_tumor_onset_female_nonsmoker,
    deterministic_tumor_onset_male_smoker,
    deterministic_tumor_onset_female_smoker
  )
)

deterministic_params[,
  `:=`(
    n_stem_cells = 10^7,
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
  )
]


use_data(deterministic_params, overwrite = TRUE)
