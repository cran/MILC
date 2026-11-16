# Archival MILC 1.0 (archival_code/) as an oracle for the validation harness.
#
# The archival functions are sourced unchanged into one environment together
# with the archival data objects. The only shim is `data()`: the archival
# functions call data(ci.lung, ..., envir = temp.MILC.Env) and then get() the
# objects from that new environment; with `data` a no-op, get() finds them by
# inheritance in the archival environment.

archival_env <- function(root = ".") {
  e <- new.env()
  for (f in list.files(file.path(root, "archival_code", "R"), full.names = TRUE)) {
    sys.source(f, envir = e)
  }
  for (f in list.files(file.path(root, "archival_code", "data"), full.names = TRUE)) {
    load(f, envir = e)
  }
  e$data <- function(..., envir) invisible(NULL)
  e
}

# Split 1:n into `n_chunks` contiguous chunks; each chunk gets its own
# L'Ecuyer-CMRG stream (parallel::nextRNGStream), so results depend on the
# seed and the number of chunks, not on the number of cores.
chunked_streams <- function(n, n_chunks, seed) {
  n_chunks <- max(1L, min(n_chunks, n))
  idx <- split(seq_len(n), cut(seq_len(n), n_chunks, labels = FALSE))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  s <- .Random.seed
  seeds <- vector("list", length(idx))
  for (k in seq_along(idx)) {
    seeds[[k]] <- s
    s <- parallel::nextRNGStream(s)
  }
  list(idx = idx, seeds = seeds)
}

par_chunks <- function(n, seed, fun, cores, n_chunks = 4L * cores) {
  cs <- chunked_streams(n, n_chunks, seed)
  res <- parallel::mclapply(seq_along(cs$idx), function(k) {
    RNGkind("L'Ecuyer-CMRG")
    assign(".Random.seed", cs$seeds[[k]], envir = globalenv())
    fun(cs$idx[[k]])
  }, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, NA, what = "try-error")
  if (any(bad)) stop("archival chunk failed: ", res[[which(bad)[1]]])
  data.table::rbindlist(res)
}

# Latent age at death from other causes: n draws per row of `strata`
# (columns smoking_status, sex, intensity_value, age).
archival_death_other <- function(e, strata, n, seed, cores) {
  grid <- strata[rep(seq_len(.N), each = n)]
  par_chunks(nrow(grid), seed, cores = cores, fun = function(ii) {
    g <- grid[ii]
    t <- vapply(seq_along(ii), function(i) {
      e$tdeath_other(
        stats::runif(1), stats::runif(1), g$smoking_status[i],
        covs_other = c(g$age[i], g$sex[i], e$d_grp(g$intensity_value[i]))
      )$time
    }, 0)
    data.table::data.table(stratum_row = g$stratum_row, time = t - g$age)
  })
}

# Latent age at death from lung cancer: n draws per row of `strata`
# (columns stage, sex, age = age at diagnosis).
archival_death_lung <- function(e, strata, n, seed, cores) {
  grid <- strata[rep(seq_len(.N), each = n)]
  par_chunks(nrow(grid), seed, cores = cores, fun = function(ii) {
    g <- grid[ii]
    t <- vapply(seq_along(ii), function(i) {
      e$tdeath_lung(
        stats::runif(1), stats::runif(1),
        c(g$stage[i], g$age[i], g$sex[i])
      )$time
    }, 0)
    data.table::data.table(stratum_row = g$stratum_row, time = t - g$age)
  })
}

# Age at onset, one draw per row of `people` (sex, ts, tq, intensity).
archival_onset <- function(e, people, seed, cores) {
  par_chunks(nrow(people), seed, cores = cores, fun = function(ii) {
    p <- people[ii]
    t <- vapply(seq_along(ii), function(i) {
      e$t_mal(stats::runif(1), p$sex[i], p$ts[i], p$tq[i], p$intensity[i])
    }, 0)
    data.table::data.table(row = ii, T_mal = t)
  })
}

# Progression from onset at age 0, n draws.
archival_progression <- function(e, par, n, seed, cores) {
  par_chunks(n, seed, cores = cores, fun = function(ii) {
    data.table::rbindlist(lapply(ii, function(i) {
      r <- e$t_prog(1, par$m, par$cdiagn, par$creg, par$cdist)
      data.table::data.table(
        T_reg = r$Treg, T_dist = r$Tdist, T_diagn = r$Tdiagn,
        D_diagn = r$Ddiagn, stage = r$stage
      )
    }))
  })
}

# End to end: archival nat_hist (`e`: archival_env() or archival_exact_env())
# for each person, redrawn until the person is alive and undiagnosed at entry
# (T_do >= entry and T_diagn >= entry). Since
# other-cause death is independent of the onset-progression-lung-death chain,
# this samples the same conditional law as the fastMILC entry condition.
# `people`: person_id, sex, status, age (entry), intensity, ts, tq.
# `par`: list by sex of list(m, cdiagn, creg, cdist).
archival_nat_hist <- function(e, people, par, seed, cores, max_tries = 1000L) {
  par_chunks(nrow(people), seed, cores = cores, fun = function(ii) {
    p <- people[ii]
    out <- vector("list", length(ii))
    for (i in seq_along(ii)) {
      q <- par[[p$sex[i]]]
      ok <- FALSE
      n_na <- 0L
      for (k in seq_len(max_tries)) {
        r <- e$nat_hist(
          c(stats::runif(5), p$age[i], p$intensity[i]), 0, p$sex[i],
          p$status[i], p$ts[i], p$tq[i], q$m, q$cdiagn, q$creg, q$cdist
        )
        # the unmodified archival lookup can return NA death times (nat_hist
        # then returns an all-NA history); such draws are redrawn and counted
        if (anyNA(c(r$T_do, r$T_dl))) {
          n_na <- n_na + 1L
          next
        }
        if (r$T_do >= p$age[i] && r$T_diagn >= p$age[i]) {
          ok <- TRUE
          break
        }
      }
      if (!ok) stop("entry condition not met for person ", p$person_id[i])
      out[[i]] <- list(
        person_id = p$person_id[i], tries = k, na_draws = n_na, T_mal = r$T_mal,
        T_reg = r$T_reg, T_dist = r$T_dist, T_diagn = r$T_diagn,
        D_diagn = r$D_diagn, stage = r$stage, T_do = r$T_do, T_dl = r$T_dl
      )
    }
    data.table::rbindlist(out)
  })
}

# Archival tables ---------------------------------------------------------------
# Knots and right-continuous CIF of an archival table element (list with
# $time, $est; duplicated rows (t, F-), (t, F+) at each step). Lung times are
# in months and are converted to years.
table_knots <- function(tab, months = FALSE) {
  t <- tab[[1]]
  F <- tab[[2]]
  tt <- sort(unique(t))
  FF <- vapply(tt, function(s) max(F[t == s]), 0)
  if (months) tt <- tt / 12
  list(t = tt, F = FF)
}

# The archival table element for a stratum, located with the archival grep.
archival_other_table <- function(e, status, agegrp, sex, intensity) {
  obj <- switch(status, never = e$never.other, former = e$former.other, e$current.other)
  key <- if (status != "current") paste(agegrp, sex) else paste(agegrp, sex, intensity)
  obj[[grep(key, names(obj))]]
}
archival_lung_table <- function(e, stage, agegrp, sex) {
  e$ci.lung[[grep(paste(stage, agegrp, sex), names(e$ci.lung))]]
}

# Exact archival law of time since entry (other causes) at times t, for entry
# age `age`: tabulated CIF (linear between knots, i.e. the archival uniform
# within a cell) plus the archival tail, residual mass 1 - F_max spread as
# Unif(0, 110 - age - h) in time since entry.
reference_cdf_other <- function(k, age, t) {
  h <- max(k$t)
  Fmax <- max(k$F)
  stats::approx(k$t, k$F, pmin(t, h), rule = 2)$y +
    (1 - Fmax) * pmin(pmax(t, 0) / (110 - age - h), 1)
}

# Exact archival law of time since diagnosis (lung cancer) at times t > 0:
# tabulated CIF with the t = 0 atom merged into the first cell (R3), linear
# between knots, residual mass 1 - F_max as a point mass at the horizon h.
reference_cdf_lung <- function(k, t) {
  h <- max(k$t)
  F <- k$F
  F[1] <- 0
  ifelse(t >= h, 1, stats::approx(k$t, F, t, rule = 2)$y)
}

# Exact-lookup archival oracle -----------------------------------------------
# A child of the archival environment in which tdeath_other and tdeath_lung
# invert the tabulated CIF exactly (linear between knots: the archival uniform
# within a cell, without the nearest-value lookup lag D8 and its NA draws).
# Everything else (nat_hist, t_mal, t_prog, grouping, grep strata, the
# other-cause tail, the residual lung point mass) is the archival code. The
# t = 0 lung atom is merged into the first cell (R3), as in fastMILC.
archival_exact_env <- function(e) {
  x <- new.env(parent = e)
  memo <- new.env()
  knots <- function(name, tab, months) {
    if (is.null(memo[[name]])) memo[[name]] <- table_knots(tab, months)
    memo[[name]]
  }
  inv <- function(k, u) {
    j <- which(k$F >= u)[1]
    if (j == 1L) return(k$t[1])
    k$t[j - 1L] + (u - k$F[j - 1L]) / (k$F[j] - k$F[j - 1L]) * (k$t[j] - k$t[j - 1L])
  }
  x$tdeath_other <- function(u1, u2, status, covs_other) {
    obj <- switch(status, never = e$never.other, former = e$former.other, e$current.other)
    age_curr <- as.numeric(covs_other[1])
    age_group <- e$age_grp(age_curr)
    key <- if (status != "current") paste(age_group, covs_other[2]) else
      paste(age_group, covs_other[2], covs_other[3])
    ind1 <- grep(key, names(obj))
    k <- knots(names(obj)[ind1], obj[[ind1]], FALSE)
    if (u1 > max(k$F)) {
      low <- age_curr + max(k$t)
      t <- u2 * (110 - low) + age_curr # archival tail
    } else {
      t <- inv(k, u1) + age_curr
    }
    list(time = t)
  }
  x$tdeath_lung <- function(u1, u2, covs_lung) {
    age_curr <- as.numeric(covs_lung[2])
    key <- paste(covs_lung[1], e$age_grp(age_curr), covs_lung[3])
    ind2 <- grep(key, names(e$ci.lung))
    k <- knots(names(e$ci.lung)[ind2], e$ci.lung[[ind2]], TRUE)
    k$F[1] <- 0 # merge the t = 0 atom into the first cell
    if (u1 > max(k$F)) {
      t <- age_curr + max(k$t) # residual point mass (archival)
    } else {
      t <- inv(k, u1) + age_curr
    }
    list(time = t)
  }
  # nat_hist must find the exact samplers: re-bind it in the child env
  x$nat_hist <- e$nat_hist
  environment(x$nat_hist) <- x
  x
}
