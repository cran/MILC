# Validation harness: fastMILC vs the archival MILC 1.0 model.
#
# Run from the package root:
#   Rscript validation/run.R [N] [cores]
# then render the report:
#   quarto render validation/report.qmd
# Writes validation/cache/results.rds (metrics and curves; no raw draws).
#
# Configuration: fastMILC runs with milc.max_age = 110 and the archival
# other-cause tail, so that every component is the archival model.
#
# References (see report.qmd for the rationale):
#   A. Mortality components: one-sample, against the exact law of each archival
#      table (tabulated CIF + archival tail; lung atom merged, residual point
#      mass). The unmodified archival samplers are compared with the same law,
#      descriptively (D8 lookup lag, NA draws).
#   B. End to end (test): fastMILC vs the archival nat_hist with an exact table
#      lookup in tdeath_other / tdeath_lung (archival_exact_env()); all other
#      archival code unchanged.
#   C. End to end (descriptive): fastMILC vs the unmodified archival nat_hist.
#   Onset and progression components: fastMILC vs archival t_mal / t_prog.

args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args) >= 1) as.numeric(args[1]) else 1e6 # end-to-end people
cores <- if (length(args) >= 2) as.integer(args[2]) else
  max(1L, min(10L, parallel::detectCores() - 1L))
n_fast_stratum <- 1e6 # fastMILC draws per mortality stratum (one-sample test)
n_arch_stratum <- 1e5 # archival sampler draws per stratum (descriptive)
n_onset <- min(N, 1e6)
n_prog <- 1e6
M <- 110
SEED <- 20261115L

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
library(data.table)
source("validation/R/archival.R")
source("validation/R/metrics.R")
options(milc.max_age = M, milc.death_other_tail = "archival")
e <- archival_env(".")
ex <- archival_exact_env(e)
timing <- list()
tic <- function() proc.time()[["elapsed"]]
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

age_mid <- c(
  "<40" = 35, "40-44" = 42.5, "45-49" = 47.5, "50-54" = 52.5, "55-59" = 57.5,
  "60-64" = 62.5, "65-69" = 67.5, "70-74" = 72.5, "75-79" = 77.5,
  "80-84" = 82.5, ">85" = 90
)
int_mid <- c("1-10" = 5, "11-20" = 15, "21+" = 30)

# A1. Other-cause death, by mortality stratum ---------------------------------
say("other-cause death components")
t0 <- tic()
so <- copy(rates_death_other$strata)[, .(stratum_row, smoking_status, agegrp, sex, intensity)]
so[, age := age_mid[agegrp]]
so[, intensity_value := fifelse(is.na(intensity), 0, int_mid[intensity])]
so[, label := paste(smoking_status, agegrp, sex, fifelse(is.na(intensity), "", intensity))]
arch_do <- archival_death_other(e, so, n_arch_stratum, SEED + 1L, cores)
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED + 2L)
comp_do <- rbindlist(lapply(seq_len(nrow(so)), function(k) {
  s <- so[k]
  kn <- table_knots(archival_other_table(e, s$smoking_status, s$agegrp, s$sex, s$intensity))
  tk <- kn$t[kn$t > 0]
  Fr <- reference_cdf_other(kn, s$age, tk)
  x <- generate_death_other(rep(s$age, n_fast_stratum), rep(s$sex, n_fast_stratum),
                            rep(s$smoking_status, n_fast_stratum),
                            rep(s$intensity_value, n_fast_stratum),
                            tail = "archival", max_age = M) - s$age
  rbind(
    knot_row("other-cause death", s$label, x, tk, Fr, "fastMILC"),
    knot_row("other-cause death", s$label, arch_do[stratum_row == s$stratum_row, time],
             tk, Fr, "archival", role = "descriptive")
  )
}))
timing$death_other <- tic() - t0

# A2. Lung-cancer death, by stratum (stage x age group at diagnosis x sex) -----
say("lung-cancer death components")
t0 <- tic()
sl <- copy(rates_death_lung$strata)[agegrp != "NA", .(stratum_row, stage, agegrp, sex)]
sl[, age := age_mid[agegrp]]
sl[, label := paste(stage, agegrp, sex)]
arch_dl <- archival_death_lung(e, sl, n_arch_stratum, SEED + 3L, cores)
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED + 4L)
comp_dl <- rbindlist(lapply(seq_len(nrow(sl)), function(k) {
  s <- sl[k]
  kn <- table_knots(archival_lung_table(e, s$stage, s$agegrp, s$sex), months = TRUE)
  h <- max(kn$t)
  # knots inside (0, h), then just before and just after the residual point
  # mass at h ((age + h) - age is not exactly h in floating point)
  tk <- c(kn$t[kn$t > 0 & kn$t < h], h * (1 - 1e-9), h * (1 + 1e-9))
  Fr <- reference_cdf_lung(kn, tk)
  x <- generate_death_cancer(rep(s$age, n_fast_stratum), rep(s$stage, n_fast_stratum),
                             rep(s$sex, n_fast_stratum)) - s$age
  rbind(
    knot_row("lung-cancer death", s$label, x, tk, Fr, "fastMILC"),
    knot_row("lung-cancer death", s$label, arch_dl[stratum_row == s$stratum_row, time],
             tk, Fr, "archival", role = "descriptive")
  )
}))
timing$death_lung <- tic() - t0

# fastMILC population ---------------------------------------------------------
say("fastMILC population, N =", N)
t0 <- tic()
options(milc.globalRNG.seed = rep(SEED, 6))
pop <- suppressWarnings(Population$new(N))$getPopulation()
timing$fast_population <- tic() - t0
f <- pop[, .(
  person_id, sex, status = smoking_status_at_spawn_age,
  intensity = smoking_intensity_at_spawn_age, entry = spawn_age,
  ts = age_smoking_start, tq = age_smoking_quit,
  T_mal = age_at_disease_onset, T_diagn = age_at_diagnosis,
  D_diagn = diameter_at_diagnosis, stage = stage_at_diagnosis,
  age_dead, cause = cause_of_death, T_dl = age_dead_from_lung_cancer
)]
f[, intgrp := intensity_group(intensity)]
people <- f[, .(person_id, sex, status, age = entry, intensity, ts, tq)]
rm(pop)

# Onset, unconditional, on the population covariates --------------------------
say("onset components, n =", n_onset)
t0 <- tic()
po <- people[seq_len(n_onset)]
arch_on <- archival_onset(e, po, SEED + 5L, cores)[order(row)]
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED + 6L)
fast_on <- generate_onset(
  params = deterministic_params[match(po$sex, deterministic_params$sex), -"sex"],
  age_smoking_start = po$ts, age_smoking_quit = po$tq,
  intensity = po$intensity, max_age = M
)
on <- data.table(
  sex = po$sex, status = po$status, intgrp = intensity_group(po$intensity),
  fast = fast_on, arch = fifelse(arch_on$T_mal < M, arch_on$T_mal, NA_real_)
)
on_groups <- unique(on[, .(sex, status, intgrp)])[order(sex, status, intgrp)]
comp_on <- rbind(
  cdf_row("onset (unconditional, censored at 110)", "all", on$fast, on$arch),
  rbindlist(lapply(seq_len(nrow(on_groups)), function(k) {
    g <- on_groups[k]
    d <- on[g, on = .(sex, status, intgrp)]
    if (nrow(d) < 1000L) return(NULL)
    cdf_row("onset (unconditional, censored at 110)",
            paste(g$sex, g$status, g$intgrp), d$fast, d$arch)
  }))
)
timing$onset <- tic() - t0

# Progression from onset ------------------------------------------------------
say("progression components, n =", n_prog)
t0 <- tic()
par_by_sex <- lapply(split(deterministic_params, by = "sex"), function(p) list(
  m = p$tumor_growth_rate,
  cdiagn = c(p$volume_diagnosis_meanlog, p$volume_diagnosis_sdlog),
  creg = c(p$volume_regional_meanlog, p$volume_regional_sdlog),
  cdist = c(p$volume_distant_meanlog, p$volume_distant_sdlog)
))
arch_pr <- archival_progression(e, par_by_sex[["male"]], n_prog, SEED + 7L, cores)
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED + 8L)
fast_pr <- generate_progression(rep(0, n_prog),
                                deterministic_params[sex == "male"][rep(1L, n_prog)])
comp_pr <- rbind(
  cdf_row("progression: time to regional disease", "onset at 0", fast_pr$age_at_regional_disease, arch_pr$T_reg),
  cdf_row("progression: time to distant disease", "onset at 0", fast_pr$age_at_distal_disease, arch_pr$T_dist),
  cdf_row("progression: time to diagnosis", "onset at 0", fast_pr$age_at_diagnosis, arch_pr$T_diagn),
  cdf_row("progression: diameter at diagnosis", "onset at 0", fast_pr$diameter_at_diagnosis, arch_pr$D_diagn),
  rbindlist(lapply(c("localized", "regional", "distant"), function(s) prop_row(
    paste("progression: stage share", s), "onset at 0",
    fast_pr$stage_at_diagnosis == s, arch_pr$stage == s
  )))
)
timing$progression <- tic() - t0

# B, C. End to end ---------------------------------------------------------------
run_arch <- function(env, seed, label) {
  say("archival nat_hist (", label, "), N =", N)
  t0 <- tic()
  a <- archival_nat_hist(env, people, par_by_sex, seed, cores)
  timing[[paste0("arch_population_", label)]] <<- tic() - t0
  a <- people[, .(person_id, sex, status, intensity, entry = age)][a, on = "person_id"]
  a[, intgrp := intensity_group(intensity)]
  acc <- a[, .(oracle = label, n = .N, mean_tries = mean(tries),
               share_redrawn_entry = mean(tries - na_draws > 1),
               na_draws_per_person = mean(na_draws)), by = .(sex, status)]
  list(a = assemble(a, M), acceptance = acc)
}
B <- run_arch(ex, SEED + 9L, "exact lookup")
C <- run_arch(e, SEED + 10L, "unmodified")

say("metrics")
e2e <- rbind(
  cbind(oracle = "exact lookup", end_to_end_metrics(f, B$a, role = "test")),
  cbind(oracle = "unmodified", end_to_end_metrics(f, C$a, role = "descriptive"))
)
curves <- rbind(
  cif_curves(f, "fastMILC"),
  cif_curves(B$a, "archival, exact lookup"),
  cif_curves(C$a, "archival, unmodified")
)

results <- list(
  meta = list(
    date = Sys.time(), N = N, n_fast_stratum = n_fast_stratum,
    n_arch_stratum = n_arch_stratum, n_onset = n_onset, n_prog = n_prog,
    M = M, seed = SEED, cores = cores,
    fastMILC_commit = tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                               error = function(e) NA_character_),
    nhppp = as.character(packageVersion("nhppp")), R = R.version.string
  ),
  mortality = rbind(comp_do, comp_dl),
  components = rbind(comp_on, comp_pr),
  end_to_end = e2e,
  curves = curves,
  acceptance = rbind(B$acceptance, C$acceptance),
  timing = timing
)
dir.create("validation/cache", showWarnings = FALSE)
saveRDS(results, "validation/cache/results.rds")
say("done")
print(rbind(
  results$mortality[, .(part = "mortality", role, verdict)],
  results$components[, .(part = "components", role, verdict)],
  results$end_to_end[, .(part = paste("e2e", oracle), role, verdict)]
)[, .N, by = .(part, role, verdict)])
