test_that("Testing make_opt_args S3 class", {
  data <- data.frame(
    treat = rep(c(1, 2, 3, 4), 25),
    y = rep(c(TRUE, FALSE), 10),
    pred = rnorm(20)
  )

  formula_test <- formula(treat ~ pred * y)

  # check data
  ## not a data.frame
  expect_error(make_opt_args(list()), "data.frame")

  # check formula
  ## no formula
  expect_error(make_opt_args(), "formula")

  ## rest is standard workflow already checked in estimate_gps()

  # check gps_methods
  expect_error(
    make_opt_args(data,
      formula_test,
      gps_method = "abc"
    ),
    "m1"
  )

  expect_error(
    make_opt_args(data, formula_test, gps_method = c("m1", "m1")),
    "Duplicates"
  )

  expect_error(make_opt_args(data, formula_test,
    gps_method = paste0("m", 1:11)
  ), "length")

  # check reference
  ## reference not in treatment levels
  expect_error(
    make_opt_args(data, formula_test, "m1", reference = "BC"),
    "unique levels"
  )


  ## reference not a string
  expect_error(
    make_opt_args(data, formula_test, "m1", reference = 3),
    "string"
  )

  # matching method
  expect_error(
    make_opt_args(data, formula_test, matching_method = "abc"),
    "matching_method"
  )

  # caliper
  ## caliper is 0
  expect_error(
    make_opt_args(data, formula_test, caliper = 0),
    "Zeros"
  )

  ## caliper is not numeric
  expect_error(
    make_opt_args(data, formula_test, caliper = "a"),
    "numeric"
  )

  ## caliper duplicated
  expect_error(
    make_opt_args(data, formula_test, caliper = c(1, 1, 1.2)),
    "Duplicates"
  )

  # order
  expect_error(
    make_opt_args(data, formula_test, order = "ab"),
    "order"
  )

  # cluster
  ## cluster exceeds number of unique treatment levels
  expect_error(make_opt_args(data, formula_test, cluster = 5), "cluster")

  # method = "nnm"
  ## validate ratio
  expect_error(
    make_opt_args(data, formula_test,
      matching_method = "nnm", ratio = 21
    ),
    "ratio"
  )

  ## validate replace
  expect_error(
    make_opt_args(data, formula_test,
      matching_method = "nnm", replace = "21"
    ),
    "replace"
  )

  ## validate ties
  expect_error(
    make_opt_args(data, formula_test,
      matching_method = "nnm", ties = "21"
    ),
    "ties"
  )

  # method = "fullopt"
  ## validate min and max controls
  expect_error(make_opt_args(data, formula_test,
    matching_method = "fullopt",
    min_controls = "21"
  ), "min_controls")

  expect_error(make_opt_args(data, formula_test,
    matching_method = "fullopt",
    max_controls = "21"
  ), "max_controls")

  # check manually number of combinations
  ## both methods
  both <- make_opt_args(data,
    formula_test,
    reference = c("1", "3"),
    gps_method = paste0("m", 1:10),
    matching_method = c("nnm", "fullopt"),
    caliper = seq(0.37, 2.21, 0.01), # 185
    order = c("desc", "random", "asc"),
    cluster = 1:3,
    replace = c(TRUE, FALSE),
    ties = FALSE,
    ratio = 1:3,
    min_controls = 1:3,
    max_controls = 3
  )

  expected_nnm <- 2 * 10 * 185 * 3 * 3 * 2 * 3
  expected_fullopt <- 2 * 10 * 185 * 3 * 3 * 3

  expect_equal(
    as.numeric(attr(both, "total_combinations")),
    expected_nnm + expected_fullopt
  )

  ## "nnm"
  nnm <- make_opt_args(data,
    formula_test,
    reference = c("1", "3", "2"),
    gps_method = c("m1", "m3", "m7", "m10"),
    matching_method = "nnm",
    caliper = seq(0.1, 2.5, 0.1), # length 25
    order = c("desc", "random"),
    cluster = 2:4,
    replace = TRUE,
    ties = FALSE,
    ratio = c(1, 2)
  )

  expected_comb <- 3 * 4 * 25 * 2 * 3 * 2

  expect_equal(as.numeric(attr(nnm, "total_combinations")), expected_comb)

  ## "fullopt"
  fullopt <- make_opt_args(data,
    formula_test,
    reference = c("1", "2", "3", "4"),
    gps_method = c("m1", "m2", "m3", "m4", "m5", "m6"),
    matching_method = "fullopt",
    caliper = seq(0.01, 3.45, 0.01), # 345
    order = "desc",
    cluster = 1:4,
    min_controls = 1:2,
    max_controls = 1:2
  )

  expected_comb <- 4 * 6 * 345 * 4 * 2 * 2

  expect_equal(as.numeric(attr(fullopt, "total_combinations")), expected_comb)
})

test_that("Optimize GPS check", {
  ## optimization code takes long to run so skip on CRAN servers
  skip_on_cran()

  # setting up data
  data <- cancer

  formula <- formula(status ~ age * sex)

  opt_args <- make_opt_args(data, formula, gps_method = "m1")

  ## test print method
  expect_no_error(print(opt_args))

  ## clean run
  suppressMessages(
    suppressWarnings(
      withr::with_seed(
        8252,
        expect_no_error(optimize_gps(data,
          formula,
          opt_args = opt_args,
          n_iter = 500
        ))
      )
    )
  )

  ## clean multicore
  suppressMessages(
    suppressWarnings(
      withr::with_seed(
        8252,
        expect_no_error(optimize_gps(data,
          formula,
          opt_args = opt_args,
          n_iter = 500,
          n_cores = 5
        ))
      )
    )
  )

  ## clean no opt_args
  suppressMessages(
    suppressWarnings(
      withr::with_seed(
        8252,
        expect_no_error(optimize_gps(data,
          formula,
          n_iter = 500
        ))
      )
    )
  )
})
