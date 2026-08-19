ARCHIVE <- "archival_code/data"
load(file.path(ARCHIVE, "current.other.rda"))
load(file.path(ARCHIVE, "former.other.rda"))
load(file.path(ARCHIVE, "never.other.rda"))

# code that parses the lists from the original MILC
parse_CIF_data <- function(
  cif,
  starting_age,
  sex_str,
  smoking_status_str,
  smoking_intensity_str,
  outcome
) {
  dat <- as.data.table(cif)
  dat[, bounds := rep(1:2, .N / 2)]
  dat[, id := floor((1:.N + 1) / 2)]
  dat <- dcast(
    dat,
    id + est ~ bounds,
    value.var = "time",
    sep = "."
  )
  new_names <- names(dat)
  new_names[which(new_names == "est")] <- "Lambda"
  new_names[which(new_names == "1")] <- "time_lb"
  new_names[which(new_names == "2")] <- "time_ub"
  setnames(dat, old = names(dat), new = new_names)
  dat[, interval := time_lb / 0.25]
  dat[, dt := time_ub - time_lb]
  dat[, dLambda := c(diff(Lambda), 0)]
  dat[, lambda := dLambda / dt]
  dat[, lambda := fifelse(lambda < 10^-8, NA, lambda)]
  dat[, (c("id", "time_ub", "time_lb", "dLambda", "dt")) := NULL]
  dat <- dat[CJ(interval = 0:40), on = "interval"]
  dat[, age := starting_age + interval * 0.25]
  dat[, lambda := nafill(lambda, type = "locf")][,
    lambda := nafill(lambda, type = "nocb")
  ]
  #dat[, Lambda := shift(cumsum(lambda * 0.25), n = 1, fill = 0, type = "lag")]
  dat[, Lambda := NULL]
  dat[, sex := sex_str]
  dat[, smoking_status := smoking_status_str]
  dat[, smoking_intensity := smoking_intensity_str]
  dat[, outcome := outcome]
  return(dat)
}

# sanitize and parse names of lists
parse_CIF_names <- function(source_list = current.other) {
  #browser()
  names_list <- names(source_list)
  fields <- names_list |>
    gsub(pattern = "dead, other", replacement = "dead_other") |>
    gsub(pattern = "cigs", replacement = "") |>
    gsub(pattern = "<40", replacement = "35-40") |>
    gsub(pattern = ">85", replacement = "85-100") |>
    gsub(pattern = "21\\+", replacement = "21-60") |>
    gsub(pattern = " +", replacement = " ") |>
    gsub(pattern = " 1$", replacement = " death_lung_cancer") |>
    strsplit(, split = " ")
  names(fields) <- names_list
  return(fields)
}

# Current smokers, other causes of death
source_list <- current.other
fields <- parse_CIF_names(source_list = current.other)
results <- list()
for (nx in names(source_list)) {
  results[[nx]] <- parse_CIF_data(
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
fields <- parse_CIF_names(source_list = former.other)
results <- list()
for (nx in names(source_list)) {
  results[[nx]] <- parse_CIF_data(
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
fields <- parse_CIF_names(source_list = never.other)
results <- list()
for (nx in names(source_list)) {
  results[[nx]] <- parse_CIF_data(
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
