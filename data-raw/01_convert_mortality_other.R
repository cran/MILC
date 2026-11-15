# Death from other causes: archival CIF tables -> nhppp step-rate inputs.
#
# Clock: years since entry. Stratum: (smoking status at entry, age group at
# entry, sex, intensity group for current smokers). The table is fixed by the
# age group at entry for the whole follow-up (archival behavior).
# Run from the package root: Rscript data-raw/01_convert_mortality_other.R

library(data.table)
source(file.path("data-raw", "00_convert_mortality_functions.R"))

ARCHIVE <- file.path("archival_code", "data")
env <- new.env()
for (f in c("never.other", "former.other", "current.other")) {
  load(file.path(ARCHIVE, paste0(f, ".rda")), envir = env)
}

# Archival names: "<agegrp> <sex> never dead, other",
# "<agegrp> <sex> former dead, other", "<agegrp> <sex> <intensity> cigs dead, other".
parse_other_names <- function(nm, smoking_status) {
  parts <- strsplit(nm, " ", fixed = TRUE)
  data.table(
    smoking_status = smoking_status,
    agegrp = vapply(parts, `[`, "", 1L),
    sex = vapply(parts, `[`, "", 2L),
    intensity = if (smoking_status == "current") {
      vapply(parts, `[`, "", 3L)
    } else {
      NA_character_
    }
  )
}

tables <- c(env$never.other, env$former.other, env$current.other)
keys <- rbindlist(list(
  parse_other_names(names(env$never.other), "never"),
  parse_other_names(names(env$former.other), "former"),
  parse_other_names(names(env$current.other), "current")
))
stopifnot(!anyDuplicated(keys))

rates_death_other <- build_rate_table(tables, keys, time_per_year = 1)

usethis::use_data(rates_death_other, overwrite = TRUE)
