#' Helper functions
#'
#' Small utility functions. Not to be exported to the user.

calendar_year_from_date <- function(date_moment) {
  as.numeric(format(date_moment, "%Y"))
}

#' @description helper function that augments
#' `test_that::expect_no_error()` to expect no error.
#' Copied from the `R6` source code.
#' @param expr Expression.
#' @return NULL
expect_no_error <- function(expr) {
  err <- FALSE
  tryCatch(force(expr), error = function(e) {
    err <<- TRUE
  })
  testthat::expect(!err, "Expected no error, but had error.")
  invisible(NULL)
}

#' Read code from text file as string
#'
#' @param codeFile Path to file
#' @return \code{codeFile} contents as a character string
read_code <- function(codeFile) {
  paste(readLines(file.path(codeFile)), collapse = "\n")
}


#' sample a category with a probability vector
#'
#' @param prob (vector or matrix) probability vector or matrix. If a matrix, probabilities are in the rows.
#' @param n number of samples to draw, if a vector. If a matrix, `n` should be the number of rows.
#' @return vector of category samples
rcat <- function(prob, n = if (is.matrix(prob)) nrow(prob) else 1) {
  if (is.vector(prob)) {
    prob <- matrix(rep(prob, each = n), nrow = n)
  }
  stopifnot(n == nrow(prob))
  n_categories <- ncol(prob)

  if (n_categories == 1) {
    return(rep(1, n))
  }

  cumsum_prop <- t(apply(prob, 1, cumsum))
  # the max standardizes the sum
  U <- runif(n = n, min = 0L, max = cumsum_prop[, n_categories])
  sampled <- rep(1, n)
  for (i in 1:(n_categories - 1)) {
    sampled <- sampled + as.numeric(U > cumsum_prop[, i])
  }
  return(sampled)
}


#' modified exponential function that is capped
#' to return a finite large value rather than `Inf`
#'
#' @param x double
#' @return double with ceiling at `getOption("kystis.exp_max")`
exp_capped <- function(x) {
  cap <- getOption("kystis.exp_max")
  exp((x < cap) * x + (x >= cap) * cap)
}

#' returns the positive part of a number or
#' `0` if the number is negative or NA
#' @param x double vector
#' @return double vector
pos_part <- function(x) {
  # attributes are copied from the first argument, so x first
  pmax(x, 0, na.rm = TRUE)
}


#' @description Ensure that the `folder` exists
#' @param folder Path to a folder. Create the path if it's not there.
#' @return Return the full path to the `folder`
ensure_folder_exists <- function(folder) {
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  return(normalizePath(folder))
}


#' @export
remove_lecuyer <- function() {
  ret <- FALSE
  if (exists(".lec.Random.seed.table", envir = .GlobalEnv)) {
    rlecuyer::.lec.CurrentStreamEnd()
    ret <- rlecuyer::.lec.exit()
  }
  invisible(ret == TRUE)
}
