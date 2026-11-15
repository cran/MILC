# Death from lung cancer: archival CIF tables -> nhppp step-rate inputs.
#
# Clock: years since diagnosis (archival tables are in months; tdeath_lung
# divides by 12). Stratum: (stage at diagnosis, age group at diagnosis, sex).
# The atom G(0) > 0 is merged into the first cell. Strata with age group "NA"
# are kept as literal labels; the archival grep lookup falls through to them
# when age_grp() returns NA.
# Run from the package root: Rscript data-raw/02_convert_mortality_lung_cancer.R

library(data.table)
source(file.path("data-raw", "00_convert_mortality_functions.R"))

ARCHIVE <- file.path("archival_code", "data")
env <- new.env()
load(file.path(ARCHIVE, "ci.lung.rda"), envir = env)

# Archival names: "<stage> <agegrp> <sex> 1".
parts <- strsplit(names(env$ci.lung), " ", fixed = TRUE)
keys <- data.table(
  stage = vapply(parts, `[`, "", 1L),
  agegrp = vapply(parts, `[`, "", 2L),
  sex = vapply(parts, `[`, "", 3L)
)
stopifnot(!anyDuplicated(keys))

rates_death_lung <- build_rate_table(env$ci.lung, keys, time_per_year = 12)

usethis::use_data(rates_death_lung, overwrite = TRUE)
