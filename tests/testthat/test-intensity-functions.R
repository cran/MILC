# Oracle values computed from the archival scalar functions
# (archival_code/R/{HT_mal,t_vol,age_grp,d_grp}.R). These pin the vectorized
# ports to the original implementation.

test_that("milc_HT_mal matches archival HT_mal (scalar cases)", {
  expect_equal(milc_HT_mal(60, "male",   20, "Y"), -0.005967192203, tolerance = 1e-9)
  expect_equal(milc_HT_mal(70, "female", NA, "N"), -0.16792642,     tolerance = 1e-8)
  expect_equal(milc_HT_mal(50, "male",   NA, "N"), -0.2183821421,   tolerance = 1e-9)
  expect_equal(milc_HT_mal(65, "female", 30, "Y"),  0.02164237074,  tolerance = 1e-9)
})

test_that("milc_HT_mal is fully vectorized (recycling)", {
  t <- c(50, 60, 70)
  out <- milc_HT_mal(t, gender = c("male", "male", "female"),
                     d = c(NA, 20, NA), smok_status = c("N", "Y", "N"))
  expect_length(out, 3)
  expect_equal(out[1], milc_HT_mal(50, "male",   NA, "N"))
  expect_equal(out[2], milc_HT_mal(60, "male",   20, "Y"))
  expect_equal(out[3], milc_HT_mal(70, "female", NA, "N"))
})

test_that("onset cumulative intensity ht_mal_int(0, t) is 0 at 0 and nondecreasing", {
  tt <- seq(0, 90, by = 1)
  Lam <- milc_ht_mal_int(0, tt, "male", 20, "Y")
  expect_equal(Lam[1], 0, tolerance = 1e-12)
  expect_true(all(diff(Lam) >= -1e-10))
})

test_that("milc_t_vol matches archival t_vol", {
  expect_equal(milc_t_vol(1e4, 4.2e-4), 9.409032888, tolerance = 1e-8)
  expect_equal(milc_t_vol(5e5, 1e-3),   6.043236876, tolerance = 1e-8)
  expect_length(milc_t_vol(c(1e4, 5e5), c(4.2e-4, 1e-3)), 2)
})

test_that("milc_age_grp matches archival age_grp", {
  ages <- c(35, 40, 42, 45, 59, 60, 85, 86, 90, 110)
  expect <- c("<40", "40-44", "40-44", "45-49", "55-59",
              "60-64", ">85", ">85", ">85", ">85")
  expect_equal(milc_age_grp(ages), expect)
  expect_true(is.na(milc_age_grp(NA)))
})

test_that("milc_d_grp matches archival d_grp", {
  ds <- c(0.5, 1, 5, 10, 11, 15, 20, 21, 40)
  expect <- c("NA", "1-10", "1-10", "1-10", "11-20", "11-20", "11-20", "21+", "21+")
  expect_equal(milc_d_grp(ds), expect)
  expect_true(is.na(milc_d_grp(NA)))
})
