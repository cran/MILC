#' Calendar year of a date
#'
#' @param date_moment `Date` vector.
#' @return numeric calendar years.
#' @keywords internal
calendar_year_from_date <- function(date_moment) {
  as.numeric(format(date_moment, "%Y"))
}

#' Remove the rlecuyer stream table
#'
#' Ends the current L'Ecuyer stream and removes `.lec.Random.seed.table` from
#' the global environment, if present.
#'
#' @return `TRUE` (invisibly) if a table was removed.
#' @export
remove_lecuyer <- function() {
  ret <- FALSE
  if (exists(".lec.Random.seed.table", envir = .GlobalEnv)) {
    rlecuyer::.lec.CurrentStreamEnd()
    ret <- rlecuyer::.lec.exit()
  }
  invisible(ret == TRUE)
}
