#' RNGStreams: class to handle RNGStreams switching for whole tasks.
#'
#' @description This class is used to handle the switching of RNGStreams for whole tasks.
#' It is based on the rlecuyer package, which provides an interface to the L'Ecuyer random number generator.
#' It generates the random number streams that are needed by Kystis (vectorized).
#' Only one instance of this class can be generated.
#' @param antithetic logical, whether to use antithetic variates.
#' @param increased_precision logical, whether to use increased precision.
#' @param stream_names character vector, names of the streams to be created. Set to package options.
#' @param seed integer vector, seed for the random number generator. Set to package options.
#' @import rlecuyer
#' @export
RNGStreams <- R6::R6Class(
  classname = "RNGStreams",
  cloneable = FALSE,
  public = list(
    initialize = function(antithetic = getOption("pneumon.antithetic"),
                          increased_precision = FALSE,
                          stream_names = getOption("pneumon.rng_stream_names"),
                          seed = getOption("pneumon.globalRNG.seed")) {
      if (exists(".lec.Random.seed.table", envir = .GlobalEnv)) {
        warning("The RNGStreams class is already initialized.
          In the old `rng` RNGStreams object, use the `rng$stopAllStreams()`
          method to remove it before creating a new object.\n
          I will erase the previous object and start a new one.")
        remove_lecuyer() # package function
      }

      stopifnot("simulation" %in% stream_names)

      private$.setGlobalSimulationSeed(seed = seed)

      private$.antithetic <- antithetic
      private$.increased_precision <- increased_precision
      private$.stream_names <- stream_names

      rlecuyer::.lec.CreateStream(stream_names)
      # .lec.Random.seed.table
      if (antithetic) {
        for (stream in stream_names) {
          rlecuyer::.lec.SetAntithetic(stream, anti = antithetic)
        }
      }
      if (increased_precision) {
        for (stream in stream_names) {
          rlecuyer::.lec.IncreasedPrecis(stream, incp = increased_precision)
        }
      }
      private$.active_stream <- "simulation"
      private$.old_kind <- rlecuyer::.lec.CurrentStream(private$.active_stream)
      invisible(self)
    },
    #' @description Print the seed table
    #' @param what character, what to print. Default is "seed_table".
    #' @return the .lec.Random.seed.table
    print = function(what = "seed_table") {
      if (what == "seed_table") {
        return(.lec.Random.seed.table)
      } else {
        stop("Unknown print option")
      }
    },
    #' @description Activate a stream. If `stream` is "R", the default RNGkind of `R` is activated.
    #' @param stream character, the name of the stream to activate, or "R" to use the default RNGkind.
    #' @return invisible(self)
    activateStream = function(stream) {
      # browser()
      if (stream == "R") {
        return(private$.unsetCurrentStream())
      }
      if (stream == private$.active_stream) {
        return(invisible(self))
      }

      stopifnot(stream %in% rlecuyer::.lec.GetStreams())
      private$.unsetCurrentStream()
      private$.setCurrentStream(stream = stream)
    },
    #' @description Get the name of the active stream.
    #' @return character, the name of the active stream. It will be "R" is no stream is active.
    getActiveStreamName = function() {
      return(private$.active_stream)
    },
    #' @description Reset the active stream.
    #' @param stream character, the name of the stream to reset. `stream` has to be an initialized stream.
    #' @return invisible(self)
    resetStream = function(stream) {
      is_active <- (stream == private$.active_stream)
      if (is_active) private$.unsetCurrentStream()
      rlecuyer::.lec.ResetStartStream(name = stream)
      if (is_active) private$.setCurrentStream(stream = stream)
      invisible(self)
    },
    #' @description Reset all streams.
    #' @return invisible(self)
    resetAllStreams = function() {
      for (stream in .lec.Random.seed.table$name) {
        self$resetStream(stream = stream)
      }
      invisible(self)
    },
    #' @description Reset the active stream to the start of the current substream.
    #' @param stream character, the name of the stream to reset. `stream` has to be an initialized stream.
    #' @return invisible(self)
    resetCurrentSubstream = function(stream) {
      is_active <- (stream == private$.active_stream)
      if (is_active) private$.unsetCurrentStream()
      rlecuyer::.lec.ResetStartSubstream(name = stream)
      if (is_active) private$.setCurrentStream(stream = stream)
      invisible(self)
    },
    #' @description Reset the active stream to the next substream.
    #' @param stream character, the name of the stream to reset. `stream` has to be an initialized stream.
    #' @return invisible(self)
    resetToNextSubstream = function(stream) {
      is_active <- (stream == private$.active_stream)
      if (is_active) private$.unsetCurrentStream()
      rlecuyer::.lec.ResetNextSubstream(name = stream)
      if (is_active) private$.setCurrentStream(stream = stream)
      invisible(self)
    },
    #' @description Stop all streams and reset the default RNGkind of `R`. (Avoid using this method.)
    #' After running this command, the object is non-functional and should be removed from the workspace.
    #' This method is useful for debugging and testing, to avoid restarting R when developing this class.
    #' @return invisible(self)
    stopAllStreams = function() {
      private$.unsetCurrentStream()
      remove_lecuyer()
      # rlecuyer::.lec.exit()
    },
    #' @description Set the object's global seed. Stops all streams,
    #' and sets the global seed for the rlecuyer package.
    #' @param seed integer vector, the seed to set.
    #' @return invisible(self)
    setSeed = function(seed) {
      self$stopAllStreams()
      self$initialize(
        antithetic = private$.antithetic,
        increased_precision = private$.increased_precision,
        stream_names = private$.stream_names,
        seed = seed
      )
      invisible(self)
    }
  ),
  private = list(
    .old_kind = NULL,
    .active_stream = NULL,
    .antithetic = NULL,
    .increased_precision = NULL,
    .stream_names = NULL,
    .setGlobalSimulationSeed = function(seed) {
      if (is.null(seed)) {
        seed <- rep(as.integer(Sys.time()), 6)
      } else if (length(seed) != 6) {
        seed <- rlecuyer::.lec.CheckSeed(seed)
      }
      rlecuyer::.lec.SetPackageSeed(seed)
      invisible(self)
    },
    .unsetCurrentStream = function() {
      private$.active_stream <- "R"
      rlecuyer::.lec.CurrentStreamEnd(kind.old = private$.old_kind)
      invisible(self)
    },
    .setCurrentStream = function(stream) {
      rlecuyer::.lec.CurrentStream(name = stream)
      private$.active_stream <- stream
      invisible(self)
    }
  )
)
