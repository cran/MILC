# Regression test for the rstream/rlecuyer symbol conflict.
#
# rstream.so and rlecuyer.so both export `user_unif_rand`; R resolves the
# symbol from the most recently loaded DLL. nhppp <= 1.1.0 imported a class
# from rstream, so `library(fastMILC)` loaded rstream after rlecuyer and the
# RNGStreams draws failed with "Stream not initialized". nhppp 1.1.0.9000 moved
# rstream to Suggests. load_all() hides the problem (load order differs), so
# the check runs the installed package in a fresh R process.

test_that("installed fastMILC simulates without loading rstream", {
  skip_on_cran()
  # skip_if_not_installed() is fooled by load_all(); look in the libraries
  installed <- file.exists(file.path(.libPaths(), "fastMILC", "DESCRIPTION"))
  skip_if(!any(installed), "fastMILC is not installed")
  rscript <- file.path(R.home("bin"), "Rscript")
  code <- paste(
    "suppressMessages(library(fastMILC))",
    "options(milc.globalRNG.seed = rep(20021115, 6))",
    "p <- Population$new(target_population_size = 50)",
    "stopifnot(nrow(p$getPopulation()) == 50)",
    "cat(\"rstream_loaded=\", \"rstream\" %in% loadedNamespaces(), \"\\n\", sep = \"\")",
    sep = "; "
  )
  out <- suppressWarnings(system2(rscript, c("-e", shQuote(code)),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  expect_true(is.null(status) || status == 0L,
              info = paste(out, collapse = "\n"))
  expect_true(any(out == "rstream_loaded=FALSE"),
              info = paste(out, collapse = "\n"))
})
