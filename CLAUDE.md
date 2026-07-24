# fastMILC — project guide & progress tracker

Refactor of the archival **MILC** (MIcrosimulation Lung Cancer) model into a
vectorized, `nhppp`-based continuous-time microsimulation, packaged for CRAN.
Also serves as a teaching example for refactoring scalar microsims to `nhppp`.

**Status:** package skeleton scaffolded and loading; mortality tables converted and
shipped as package data; ported kernels tested (25/25 passing). Event generators are
documented stubs pending D1/D2 implementation. Open items: D1 upper tail, new D6
(t=0 atom for lung death), F=1 infinite-rate handling.

---

## 1. Provenance

- `archival_code/` — the original CRAN package **MILC 1.0** (Chrysanthopoulou, 2014).
  Per-individual, scalar sampling. This is the *reference implementation*; the
  refactor must reproduce its distributions (see §6 validation).
- `example_to_imitate/` — a generic-cancer teaching model (`ToyModel` R6 class)
  that already uses `nhppp`. This is the *architectural template* we imitate:
  population `data.table` + R6 encapsulation + vectorized `nhppp::vdraw_*`.
- Environment: R 4.5.3, `nhppp` 1.0.5, `data.table`, `R6`.
- `nhppp` is **authored by the maintainer (TAT)**; source at
  https://github.com/cran/nhppp.git — treat its vignettes/docs as the
  authoritative reference for sampler semantics (`vdraw_*`, majorizers, step fns).

## 2. The MILC model, event by event

`nat_hist()` stitches four stochastic components per individual:

| Component (old fn)   | Object                                            | Refactor target |
|----------------------|---------------------------------------------------|-----------------|
| `t_mal`              | First arrival of NHPPP; **closed-form Λ** = `HT_mal` (TSCE-type), piecewise over smoking segments never→`ts`→`tq` | `nhppp::vdraw_cumulative_intensity` (invert Λ; NOT majorizer) |
| `t_prog`             | Deterministic Gompertz tumor growth (`t_vol`) + random lognormal volume thresholds at dx/regional/distant, with rejection for ordering & bounds | Vectorized rejection/transform sampler — **not** a Poisson process |
| `tdeath_other`       | Inverse-CDF from tabulated **CIF** step functions, strata = age-group × sex (× smoking intensity for current smokers) | See decision D1 |
| `tdeath_lung`        | Inverse-CDF from tabulated CIF `ci.lung`, strata = stage × age-group × sex | See decision D1 |

Downstream logic in `nat_hist`: min(`T_do`, `T_dl`) → cause & age at death;
prediction horizon `T_pred`; incidence/death indicators; exclusion of cases with
`T_death < T_curr`.

### Data objects (archival, in `archival_code/data/*.rda`)
- `ci.lung` — list of 68; keys `"<stage> <agegrp> <sex> 1"`; each `$time,$est,$var` (CIF).
- `never.other` / `former.other` / `current.other` — lists (22 / … ); keys
  `"<agegrp> <sex> [intensity] <status>"`; each `$time,$est,$var` (CIF).
- Smoking intensity groups: `d_grp()` → NA / 1-10 / 11-20 / 21+.
- Age groups: `age_grp()` → `<40 … >85` in 5-yr bands.

## 3. Target architecture (imitating `ToyModel`)

- `data.table` population, one row per individual, all params as columns.
- R6 class `MILC` (name TBD, D4): `initialize / update / reset / getPopulation`,
  private per-event generators, `copy = TRUE` semantics, `setindex` on join keys.
- Pure vectorized generator functions (one per event), each thin wrappers over
  `nhppp::vdraw_*`, taking a population `data.table` + rate/param inputs.
- Parameters and reference tables shipped as package data (`data/`), regenerated
  by scripts in `data-raw/` (not shipped).

## 4. Proposed package layout (CRAN-oriented)

```
fastMILC/
├── DESCRIPTION            # Imports: nhppp, data.table, R6
├── NAMESPACE              # roxygen2-generated (no exportPattern)
├── R/
│   ├── milc-class.R       # R6 MILC class
│   ├── intensity-fns.R    # HT_mal / Lambda + inverse, vectorized
│   ├── generate-onset.R   # t_mal  -> vdraw_cumulative_intensity
│   ├── generate-progression.R
│   ├── generate-death-other.R
│   ├── generate-death-cancer.R
│   └── milc-package.R      # @keywords internal doc, imports
├── data/                  # *.rda reference tables (CIF or rate matrices per D1)
├── data-raw/              # scripts that build data/
├── tests/testthat/        # unit + distributional-equivalence tests
├── vignettes/             # teaching vignette: scalar -> nhppp refactor
└── man/                   # roxygen2-generated
```

## 5. DECISIONS (resolved 2026-07-13)

- **D1 — mortality tables. RESOLVED.** Do NOT keep a bespoke inverse-CDF sampler.
  Re-express the other-cause and lung-cancer death tables as **regular step-intensity
  functions** and sample with `vdraw_sc_step_regular`; where a stratum's table has
  only a few steps, instead hand-write a closed-form `Lambda`/`Lambda_inv` for it and
  use `vdraw_cumulative_intensity`. nhppp-native throughout; validate induced
  marginals against the archival tables. (Note: original CIF tables are competing-risks
  cumulative incidence on irregular time — the conversion to cause-specific hazards on
  a regular grid is the modeling step to get right and document.)
  **Competing-risks semantics: match the archival model exactly** — independent
  latent event times ("time bombs"), one per cause, sampled from each cause's own
  tabulated (marginal/net) curve, then take the minimum. Do NOT "upgrade" to genuine
  cause-specific sub-distribution hazards during the port; that is a behavior change,
  not a refactor. The step-intensity tables therefore represent each cause's marginal
  latent-time distribution. This preserves old↔new equivalence for validation (§6).

  **Tabulated-part conversion (verified on `never.other` `<40 female`):** collapse the
  duplicated (t,F⁻)/(t,F⁺) step rows to a right-continuous CIF `F(t)` at unique
  breakpoints (0.25-yr grid, time-since-entry); latent survival `S=1-F`; step hazard
  per cell `λ_k = (log S(t_k) - log S(t_{k+1}))/Δ_k`. First-arrival of the NHPPP then
  reproduces `F` at every breakpoint to fp precision. Clock is time-since-entry keyed
  by **entry** age-group×sex(×intensity); table does NOT switch as the person ages —
  preserve this.

  **Tail decision (residual ~99% with no event in the tabulated window): "match
  distribution, fix quirks".** The archival tail smears residual mass as
  `Uniform(age_curr, 94)` via `u2*(up-low)+age_curr` — two apparent bugs: it adds
  `age_curr` instead of `low` (overlaps the tabulated window) and the upper bound
  collapses to 94 not 110. Refactor uses the corrected **`Uniform(age_curr+16, 110)`**
  (no overlap, bound at 110). Faithful in spirit, NOT bit-identical → validation
  tolerances loosen in the tail only. A `Uniform(a,b)` residual has increasing hazard
  `1/(b-t)`, so tail cells are per-individual (bound depends on entry age); note the
  16 is the `<40`-group horizon and is group-specific (`max(time)` per stratum).
- **D2 — onset sampler. RESOLVED.** `vdraw_cumulative_intensity`, inverting the
  closed-form Λ = `HT_mal` (piecewise over smoking segments). No majorizer.
- **D3 — scope. RESOLVED.** MILC (lung cancer) only. **No generalization** to a
  generic natural-history skeleton. Keep the code concrete and MILC-specific; the
  reusable abstraction lives in the `nhppp-refactor` skill (§7), not in the package.
- **D4 — names. RESOLVED.** Package `fastMILC`, R6 class `MILC`, generator functions
  `generate_*` as in the layout below. Preserve MILC terminology.
- **D5 — reproducibility. RESOLVED.** Use **RngStreams** (L'Ecuyer) for seeding,
  via `nhppp::rng_stream_*` wrappers. Consistent with the archival code, which already
  carried L'Ecuyer stream state (`.lec.Random.seed.table` in `MILC-internal.R`). One
  named stream per stochastic component so old↔new comparisons and teaching runs are
  reproducible and parallel-safe.

- **D6 — t=0 point mass in lung-death tables. OPEN (new, found during scaffold).**
  The `ci.lung` CIFs start at `F(0) > 0` (distant up to 0.137, regional 0.003–0.056,
  localized up to 0.028) — an instantaneous probability of death at diagnosis. All
  other-cause tables have `F(0)=0`. The shipped step-hazard cells therefore encode the
  distribution **conditional on surviving the atom**; `F0` is carried per stratum in
  `lung_strata`/`other_strata`. `generate_death_cancer` must sample the atom as a
  Bernoulli(`F0`) immediate event, then the step process conditional on survival.
  Analogous to the archival inverse-CDF, which interpolated over the `F(0)` value
  directly. Decide the exact split when implementing the lung-death generator.

## 6. Validation strategy (scientific correctness ≠ software correctness)

- Per-stratum KS tests: refactored vs. archival marginals for `T_mal`, `T_do`,
  `T_dl`, stage distribution, diameter at diagnosis.
- Aggregate checks: lung-cancer incidence, cause-specific mortality, exclusion rate.
- Seed-controlled reproducibility of the vectorized path.
- Keep archival scalar functions available in `tests/` as an oracle.

## 7. Abstractable skills (secondary deliverable)

Distill reusable Claude Code skills from this refactor:
1. **`nhppp-refactor`** *(primary)* — classify a scalar event-time sampler →
   {closed-form Λ → `vdraw_cumulative_intensity`; only λ(t) → `vdraw_intensity` +
   `get_step_majorizer`; tabulated/step rate → `vdraw_sc_step_regular`} → move
   per-individual params to a `data.table` of `lambda_args` → replace
   `uniroot`/rejection with one vectorized draw → validate equivalence. MILC is its
   worked example.
2. **`r6-microsim-scaffold`** — the `ToyModel` idiom as a template.
3. **`microsim-validation`** — old-vs-new distributional equivalence harness.

## 8. Progress log

- 2026-07-13 — Read archival code + nhppp example; drafted architecture, mapping,
  and this tracker. Resolved D1–D5: MILC-specific only; nhppp-native mortality via
  regular step-intensity (or closed-form Λ for few-step strata); onset via
  `vdraw_cumulative_intensity`; RngStreams for seeding; names `fastMILC`/`MILC`.
  No code written yet.
- 2026-07-13 — Converted all four mortality tables to tidy step-hazard data.tables
  (§8b); scaffolded the `fastMILC` package (§8c). `R CMD check` clean (0E/0W/3 trivial
  notes). Ported kernels validated against archival scalars; 25/25 tests pass. Found
  the t=0 lung-death atom → new decision D6. Generators remain stubs pending D1/D2/D6.

## 8b. Converted mortality data.tables (2026-07-13)

`data-raw/01_convert_mortality_tables.R` converts all four archival CIF list objects
to tidy step-hazard `data.table`s over the **tabulated horizon only** (upper tail
deferred — TAT decides separately). Outputs in `data-raw/converted/*.rds`:
- `other_cells` (4074 rows): `status, agegrp, sex, intensity, t_lo, t_hi, F_lo, F_hi, lambda`
  — combines never/former (intensity = NA) + current (intensity ∈ 1-10/11-20/21+).
- `lung_cells` (6617 rows): `stage, agegrp, sex, t_lo, t_hi, F_lo, F_hi, lambda`;
  time converted months→years (÷12).
- `other_strata`, `lung_strata`: per-stratum `horizon_yr, F_max, residual_mass, n_cells, any_inf_rate`.

Key facts / edges surfaced (all fold into the pending tail decision):
- Horizons are **stratum-specific**: never/former = 16 yr; current mostly 10 yr (a few
  3.75–9.25); lung ≈ 20.9 yr.
- **F=1 strata → Inf terminal hazard** (certain event at horizon, *opposite* of residual
  mass): `current >85 female 21+` (3.75 yr) and the fallback `distant NA male`. Need a
  finite-rate representation (cap, or explicit point mass) — decide with the tail.
- **F=0 degenerate:** `localized NA female` (all-zero CIF → latent time +Inf → 0 cells).
- **NA-agegrp fallback strata** (`distant NA male`, `localized NA female`) exist because
  the archival `grep` matches them only when computed `age_grp(age)` is NA. Preserved
  as literal `agegrp = "NA"`.

## 8c. Package scaffold (2026-07-13)

Skeleton built at repo root; `R CMD check` = 0 errors / 0 warnings / 3 notes (notes:
build-time timestamp check [env], `nhppp` unused [stubs not yet wired], `.claude`
hidden dir [now .Rbuildignore'd]). Structure:
- `R/fastMILC-package.R` — package doc, imports, globalVariables.
- `R/intensity-functions.R` — ported vectorized kernels `milc_HT_mal`, `milc_ht_mal_int`,
  `milc_t_vol`, `milc_age_grp`, `milc_d_grp` (validated vs archival scalars).
- `R/generators.R` — `generate_{onset,progression,death_other,death_cancer}` STUBS
  (`stop()` naming the blocking decision).
- `R/milc-class.R` — R6 `MILC` (initialize/update/reset/getPopulation) + private
  generator dispatch + `.resolve_trajectory` stub.
- `R/data.R` — dataset docs. `data/` holds `other_cells/other_strata/lung_cells/lung_strata`.
- `data-raw/01_convert_mortality_tables.R` (converter) + `02_make_package_data.R` (→ data/).
- `tests/testthat/` — `test-intensity-functions.R` (oracle vs archival),
  `test-mortality-data.R` (CIF round-trip incl. F0 atom; hazard well-formedness). 25/25 pass.
- `vignettes/fastMILC.Rmd` — placeholder.

**Immediate next steps:** implement `generate_death_other`/`generate_death_cancer`
(needs D1 tail + D6 atom + F=1 rate handling), then `generate_onset` (D2), then
`generate_progression`, then `.resolve_trajectory` (nat_hist stitching). Wiring nhppp
into the generators clears the unused-Imports note.

## 9. Working agreements

- Plan first; confirm D1–D5 before writing package code.
- Write tests (happy path + edge) for each generator; run `R CMD check` /
  `testthat` before any commit. Don't commit unless asked.
- Preserve MILC terminology and the maintainer's naming where technically defensible.
- Distinguish exact reproduction from approximation explicitly in code comments and
  the vignette.
