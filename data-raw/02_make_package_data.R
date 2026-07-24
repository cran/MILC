# data-raw/02_make_package_data.R
# Promote the converted mortality tables to package data/ objects.
# Run from the package root:  Rscript data-raw/02_make_package_data.R
# Regenerate whenever 01_convert_mortality_tables.R or the archival tables change.

source("data-raw/01_convert_mortality_tables.R")   # defines other_cells, other_strata, lung_cells, lung_strata

usethis::use_data(
  other_cells, other_strata, lung_cells, lung_strata,
  overwrite = TRUE, version = 3
)
