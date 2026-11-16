# fastMILC

`fastMILC` is a vectorized re-implementation of MILC, the continuous-time
microsimulation model of the natural history of lung cancer
(Chrysanthopoulou, 2014; *International Journal of Microsimulation*, 2017).
For each person it simulates the onset of the first malignant cell (two-stage
clonal expansion over the smoking history), Gompertz tumor growth, the ages at
regional and distant disease, the age, stage and tumor diameter at clinical
diagnosis, and death from lung cancer or from other causes.

The archival MILC package draws each person's history with scalar inverse-CDF
and root-finding steps. `fastMILC` draws each component for the whole
population at once with the [`nhppp`](https://github.com/bladder-ca/nhppp)
samplers over a `data.table`, and uses named L'Ecuyer random number streams,
one per model component.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ttrikalin/fastMILC", ref = "refactor")
```

`fastMILC` requires `nhppp` (>= 1.1.1), which `install_github()` installs
from [bladder-ca/nhppp](https://github.com/bladder-ca/nhppp).

## Example

```r
library(fastMILC)
library(data.table)

rng <- RNGStreams$new(seed = rep(20021115, 6))

# 5,000 people drawn from the synthetic NLST-like recipe, simulated from
# study entry with the baseline parameters
pop <- Population$new(target_population_size = 5000, rng_stream = rng)$getPopulation()

# lung-cancer diagnoses and deaths over the lifetime
pop[, .(diagnosed = mean(!is.na(age_at_diagnosis)),
        lung_cancer_deaths = mean(cause_of_death %in% "lung_cancer"))]

# stage at diagnosis
pop[!is.na(stage_at_diagnosis), .N, by = stage_at_diagnosis]

# a scenario: 20% faster tumor growth, same random number streams
fast <- copy(deterministic_params)[, tumor_growth_rate := 1.2 * tumor_growth_rate]
pop_fast <- Population$new(5000, parameters = fast, rng_stream = rng)$getPopulation()
```

The vignette shows more: the simulated columns, incidence within windows after
entry, reproducibility, scenarios, custom recipes, the package options and the
component functions:

```r
vignette("fastMILC", package = "fastMILC")
```

(source: [`vignettes/fastMILC.Rmd`](vignettes/fastMILC.Rmd)).

## Model and differences from the archival code

The model mathematics and every difference from the archival code are
documented in [`model_description/model.tex`](model_description/model.tex).
The main differences are intended:

- **Maximum age.** 100 by default (`options(milc.max_age = ...)`); the archival
  code used 110 for other-cause death and 150 for onset. With the archival
  other-cause tail, $\mathrm{Unif}(\text{entry},\ M - h)$, a lower maximum age
  moves other-cause deaths earlier: for the NLST-like recipe, 25% die of other
  causes within 5 years of entry at $M = 100$ against 19% at $M = 110$.
- **Entry condition.** Only people alive and without a diagnosed lung cancer at
  study entry enter the simulation; the archival code flagged the others for
  exclusion.
- **Exact sampling of the mortality tables.** The archival lookup lags its own
  tables (see below); `fastMILC` reproduces them at every knot.
- **Population sampling.** People are drawn from the recipe with probability
  proportional to the recipe weight.

## Validation against the archival code

The harness in [`validation/`](validation/) compares `fastMILC` with the
archival MILC 1.0 package (`Rscript validation/run.R`, then
`quarto render validation/report.qmd`). For the comparison `fastMILC` runs with
maximum age 110 and the archival other-cause tail, so that every component
specifies the archival model. Both sides simulate the same 1,000,000 people
from the recipe, conditional on being alive and undiagnosed at entry.

**Criterion.** Effect sizes with tolerance 0.005: the sup-norm distance between
distribution functions (cumulative incidence functions where events compete)
and the absolute difference of proportions. Each estimate has an approximate
95% Monte Carlo bound $b$. A comparison is *equivalent* if
$|\hat\delta| + b \le 0.005$, *different* if $|\hat\delta| - b > 0.005$, and
*inconclusive* otherwise.

| Comparison | Reference | Result |
|---|---|---|
| Mortality tables, 176 strata (110 other-cause, 66 lung) | exact law of each archival table | all equivalent; largest distance 0.0014 |
| Onset of the first malignant cell | archival `t_mal` | equivalent overall (0.0006); 8 subgroups inconclusive, distances ≤ 0.003 |
| Progression: times to regional and distant disease and to diagnosis, diameter, stage | archival `t_prog` | all 7 equivalent; distances ≤ 0.0012 |
| End to end, whole population and 9 subgroups (198 comparisons) | archival `nat_hist` with an exact table lookup | 111 equivalent, 87 inconclusive, 0 different |

For the whole population, the cumulative incidence of diagnosis, of
lung-cancer death and of other-cause death against time since entry, and the
proportions diagnosed and dead of each cause within 5, 10 and 20 years and over
the lifetime, are all equivalent; for example, 4.78% (`fastMILC`) and 4.75%
(archival) are diagnosed, and 3.58% and 3.54% die of lung cancer. The
inconclusive comparisons are subgroups and quantities conditional on diagnosis
(about 48,000 diagnoses per side), where the Monte Carlo bound exceeds the
tolerance.

**Findings about the archival code.** The archival mortality samplers do not
reproduce their own tables, which is why the end-to-end test uses an archival
`nat_hist` in which only the table lookup is exact:

- The lookup selects a cell by its position in the de-duplicated vector of
  distances `unique(abs(u1 - est))` and uses that position in `unique(time)`.
  The tables store each step twice, so the draws land one to two cells late:
  0.002--0.006 (median) from the table for other-cause death, and up to 0.27
  for death after a distant-stage diagnosis.
- When the nearest tabulated value is the last one, the lookup returns `NA`
  (up to 0.1% of draws in 50 of the 110 other-cause tables, 25% in one); the
  archival `nat_hist` then returns an all-`NA` history.

Against the unmodified archival `nat_hist`, population incidence and mortality
differ by at most 0.003, and survival after diagnosis by 0.03 (localized), 0.07
(regional) and 0.19 (distant stage).

## References

- S. A. Chrysanthopoulou. *MILC: MIcrosimulation Lung Cancer model.* R package
  version 1.0, 2014.
- S. A. Chrysanthopoulou. MILC: A microsimulation model of the natural history
  of lung cancer. *International Journal of Microsimulation*, 10(3):5--26, 2017.

## License

GPL (>= 2).
