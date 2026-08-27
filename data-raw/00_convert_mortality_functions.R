# code that parses the lists from the original MILC
parse_CIF_data_other_deaths <- function(
  cif,
  starting_age,
  sex_str,
  smoking_status_str,
  smoking_intensity_str,
  outcome,
  step_size = 0.25
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
  dat[, interval := time_lb / step_size]
  dat[, dt := time_ub - time_lb]
  dat[, dLambda := c(diff(Lambda), 0)]
  dat[, lambda := dLambda / dt]
  dat[, lambda := fifelse(lambda < 10^-8, NA, lambda)]
  dat[, (c("id", "time_ub", "time_lb", "dLambda", "dt")) := NULL]
  dat <- dat[CJ(interval = 0:40), on = "interval"]
  dat[, age := starting_age + interval * step_size]
  dat[, lambda := nafill(lambda, type = "locf")][,
    lambda := nafill(lambda, type = "nocb")
  ]
  #dat[, Lambda := shift(cumsum(lambda * step_size), n = 1, fill = 0, type = "lag")]
  dat[, Lambda := NULL]
  dat[, sex := sex_str]
  dat[, smoking_status := smoking_status_str]
  dat[, smoking_intensity := smoking_intensity_str]
  dat[, outcome := outcome]
  return(dat)
}

# sanitize and parse names of lists
parse_CIF_names_other_deaths <- function(source_list = current.other) {
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


parse_CIF_names_lung_cancer_deaths <- parse_CIF_names_other_deaths

parse_CIF_data_lung_cancer_deaths <- function(
  cif,
  starting_age,
  sex_str,
  smoking_status_str,
  smoking_intensity_str,
  outcome
) {
  ### time appears to be in days?
  ### there appears to be a 0-0 day mortality (operative?);
  ### we will add that Lambda_{0,0} to the Lambda_{0,1}
  cif$time <- cif$time[3:length(cif$time)]
  cif$est <- c(rep(sum(cif$est[1:3]), 2), cif$est[5:length(cif$est)])
  cif$var <- NULL
  dat <- parse_CIF_data_other_deaths(
    cif = cif,
    starting_age = starting_age,
    sex_str = sex_str,
    smoking_status_str = smoking_status_str,
    smoking_intensity_str = smoking_intensity_str,
    outcome = outcome,
    step_size = 1 / 365
  )
  return(dat)
}
