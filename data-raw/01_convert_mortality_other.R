ARCHIVE <- "archival_code/data"
load(file.path(ARCHIVE, "current.other.rda"))
load(file.path(ARCHIVE, "former.other.rda"))
load(file.path(ARCHIVE, "never.other.rda"))

source(file.path('data-raw', "00_convert_mortality_functions.R"))

# Current smokers, other causes of death
source_list <- current.other
fields <- parse_CIF_names_other_deaths(source_list = source_list)
results <- list()
for (nx in names(source_list)) {
  results[[nx]] <- parse_CIF_dat_other_deaths(
    cif = source_list[[nx]],
    starting_age = as.numeric(substr(fields[[nx]][1], 1, 2)),
    sex_str = fields[[nx]][2],
    smoking_status_str = "current",
    smoking_intensity_str = fields[[nx]][3],
    outcome = "death_from_other_causes"
  )
}
OC <- rbindlist(results)


# Former smokers, other causes of death
source_list <- former.other
fields <- parse_CIF_names_other_deaths(source_list = source_list)
results <- list()
for (nx in names(source_list)) {
  results[[nx]] <- parse_CIF_data_other_deaths(
    cif = source_list[[nx]],
    starting_age = as.numeric(substr(fields[[nx]][1], 1, 2)),
    sex_str = fields[[nx]][2],
    smoking_status_str = "former",
    smoking_intensity_str = fields[[nx]][3],
    outcome = "death_from_other_causes"
  )
}
OF <- rbindlist(results)

## never smokers, other causes of death
source_list <- never.other
fields <- parse_CIF_names_other_deaths(source_list = source_list)
results <- list()
for (nx in names(source_list)) {
  results[[nx]] <- parse_CIF_data_other_deaths(
    cif = source_list[[nx]],
    starting_age = as.numeric(substr(fields[[nx]][1], 1, 2)),
    sex_str = fields[[nx]][2],
    smoking_status_str = "never",
    smoking_intensity_str = fields[[nx]][3],
    outcome = "death_from_other_causes"
  )
}
ON <- rbindlist(results)

#### Other causes of death:
#### The original MILC notation conflates age and time intervals.
#### Thus the tables include different information for the same age
#### when all are put together.
#### We will use all overlapping information on ages; and average
#### the entries where available
death_other <- rbindlist(list(OF, OC, ON))
death_other[, interval := NULL]
setkey(death_other, age, sex, smoking_status, smoking_intensity)
death_other[,
  group_id := .GRP,
  by = .(age, sex, smoking_status, smoking_intensity)
]
death_other <- death_other[,
  .(
    age = first(age),
    sex = first(sex),
    smoking_status = first(smoking_status),
    smoking_intensity = first(smoking_intensity),
    outcome = first(outcome),
    lambda = mean(lambda, rm.na = TRUE)
  ),
  by = group_id
]

# death from other causes by 0.25 years of age
death_other_x0.25 <- dcast(
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
