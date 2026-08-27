#### THIS IS WRONG
####  -- JUST FIX FROM TIME OF DIAGNOSIS? PROBABLY
####  --- DISEASE ONSET? NOT LIKELY -- SEE OLD CODE

ARCHIVE <- "archival_code/data"
load(file.path(ARCHIVE, "ci.lung.rda"))

source(file.path('data-raw', "00_convert_mortality_functions.R"))

source_list <- ci.lung
fields <- parse_CIF_names_lung_cancer_deaths(source_list = source_list)


results <- list()
for (nx in names(source_list)) {
  results[[nx]] <- parse_CIF_data_lung_cancer_deaths(
    cif = source_list[[nx]],
    starting_age = as.numeric(substr(fields[[nx]][2], 1, 2)),
    sex_str = fields[[nx]][3],
    smoking_status_str = "any",
    smoking_intensity_str = "any",
    outcome = "death_from_lung_cancer"
  )
  results[[nx]][, stage := fields[[nx]][1]]
}
death_LC <- rbindlist(results)

### Warnings:
## Some strata are not described -- see for example:
# > source_list[["distant NA male 1"]]
# $time
# [1] 0 5 5 5

# $est
# [1] 0 0 1 1

# $var
# [1] 0 0 1 1

drop_vars <- c("interval", "smoking_status", "smoking_intensity")
death_LC[, (drop_vars) := NULL]
death_LC <- death_LC[!is.na(age) | !is.na(lambda)]
setkey(death_LC, age, sex, stage)
death_LC[,
  group_id := .GRP,
  by = .(age, sex, stage)
]
death_LC <- death_LC[,
  .(
    age = first(age),
    sex = first(sex),
    stage = first(stage),
    outcome = first(outcome),
    lambda = mean(lambda, rm.na = TRUE)
  ),
  by = group_id
]

# death by day
death_LC_xday <- dcast(
  death_other,
  sex + smoking_status + smoking_intensity ~ age,
  value.var = 'lambda'
)

# death from other causes by 1 years of age
death_other[, age := floor(age)][, group_id := NULL]
setkey(death_other, age, sex, smoking_status, smoking_intensity)


death_other <- death_other[,
  .(
    # age = first(age),
    # sex = first(sex),
    # smoking_status = first(smoking_status),
    # smoking_intensity = first(smoking_intensity),
    outcome = first(outcome),
    lambda = mean(lambda, rm.na = TRUE)
  ),
  by = .(age, sex, smoking_status, smoking_intensity)
]

death_other_x1 <- dcast(
  death_other,
  sex + smoking_status + smoking_intensity ~ age,
  value.var = 'lambda'
)

use_data(death_other_x1, death_other_x0.25, overwrite = TRUE)
