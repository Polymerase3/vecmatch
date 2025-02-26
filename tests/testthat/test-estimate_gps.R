## --testing formals: formula & data -------------------------------------------
test_that("Formals checking: formula", {
  data <- data.frame(
    y = runif(20),
    treat2 = c(runif(19), NA),
    group = rep(c(TRUE, FALSE), 10),
    sex = rep(c("M", "F", "F", "M"), 5)
  )
  treat <- rep(c(0, 1), each = 10)
  treat_fail <- c(rep(c(0, 1), 9), c(NA, NA))
  pred <- runif(20)
  pred_fail <- runif(21)

  expect_error(estimate_gps(), regexp = "missing")
  expect_error(estimate_gps(NULL), regexp = "treatment")
  expect_error(estimate_gps(treat ~ pred_fail), regexp = "samples")
  expect_error(estimate_gps(treat_fail ~ pred), regexp = "NA")
  expect_warning(estimate_gps(data$y ~ data$group))
  expect_no_error(estimate_gps(treat ~ pred))
})

## --testing formals: method, ref, logicals-------------------------------------
test_that("Formals checking:  method, ref, logicals", {
  data <- data.frame(
    treat = rep(c(1, 2, 3, 4, 5), 20),
    y = rep(c(TRUE, FALSE), 10),
    pred = runif(20)
  )

  # method
  expect_error(estimate_gps(y ~ pred, data, method = c()), regexp = "string")
  expect_error(estimate_gps(y ~ pred, data, method = "error"),
    regexp = "method"
  )
  expect_no_error(estimate_gps(y ~ pred, data, method = NULL))
  expect_no_error(estimate_gps(y ~ pred, data, method = "vglm"))

  # ref
  expect_error(estimate_gps(y ~ pred, data, method = NULL, reference = TRUE),
    regexp = "reference"
  )
  expect_error(estimate_gps(y ~ pred, data, method = NULL, reference = "FAIL"),
    regexp = "unique"
  )
  expect_warning(estimate_gps(treat ~ pred, data,
    method = "multinom",
    reference = "1",
    ordinal_treat = c(1, 3, 2, 5, 4)
  ))
  expect_no_error(estimate_gps(y ~ pred, data, method = NULL, reference = NULL))
  expect_no_error(estimate_gps(y ~ pred, data,
    method = "multinom",
    reference = "TRUE"
  ))

  # logicals
  expect_error(estimate_gps(y ~ pred, data, fit_object = "fail"),
    regexp = "flag"
  )

  expect_output(estimate_gps(y ~ pred, data, verbose_output = TRUE))
  expect_no_error(estimate_gps(y ~ pred, data, fit_object = TRUE))
})

## --testing formals: missing, by-----------------------------------------------
test_that("Formals checking: missing and by", {
  data <- data.frame(
    treat = rep(c("A", "B", "C"), 7),
    y = runif(21),
    group = rep(c(TRUE, FALSE, TRUE), 7),
    sex = c(rep(c("M", "F", "F", "M"), 5), "M")
  )

  expect_no_error(estimate_gps(treat ~ y, data,
    method = "multinom",
    by = "sex"
  ))
  expect_no_error(estimate_gps(treat ~ y, data,
    method = "vglm",
    by = "sex"
  ))
  expect_error(estimate_gps(treat ~ y, data, method = "vglm", by = "abc"),
    regexp = "stratify"
  )
})

## --testing formals: link------------------------------------------------------
test_that("Formals checking: link", {
  data <- data.frame(
    treat = rep(c("A", "B", "C"), 7),
    pred = runif(21)
  )

  expect_error(
    estimate_gps(treat ~ pred, data,
      method = "multinom",
      link = list()
    ),
    regexp = "string"
  )

  expect_error(
    estimate_gps(treat ~ pred, data,
      method = "multinom",
      link = "fail_link"
    ),
    regexp = "link"
  )

  expect_no_error(estimate_gps(treat ~ pred, data,
    method = "multinom",
    link = "generalized_logit"
  ))
})

## --testing formals: ordinal_treat---------------------------------------------
test_that("Formals checking: ordinal_treat", {
  data <- data.frame(
    treat = rep(c(1, 2, 3, 4, 5), 20),
    pred = runif(100)
  )
  data$treat <- factor(data$treat, levels = c(1, 2, 3, 4, 5), ordered = TRUE)

  expect_error(
    estimate_gps(treat ~ pred, data,
      method = "multinom",
      ordinal_treat = list(1)
    ),
    regexp = "atomic"
  )

  expect_error(
    estimate_gps(treat ~ pred, data,
      method = "multinom",
      ordinal_treat = c(1, 2, 3)
    ),
    regexp = "levels"
  )

  expect_no_error(estimate_gps(treat ~ pred, data,
    method = "multinom",
    ordinal_treat = c(1, 3, 2, 5, 4)
  ))
})

## --testing formals: subset----------------------------------------------------
test_that("Formals checking: subset", {
  data <- data.frame(
    treat = rep(c(1, 2, 3, 4, 5), 20),
    pred = runif(100),
    subset_pass = rep(c(TRUE, FALSE), 50),
    subset_fail = rep(10, 100)
  )

  expect_error(
    estimate_gps(treat ~ pred, data,
      method = "multinom",
      subset = list()
    ),
    regexp = "string"
  )

  expect_error(
    estimate_gps(treat ~ pred, data,
      method = "multinom",
      subset = "some_col"
    ),
    regexp = "provided dataset"
  )

  expect_error(
    estimate_gps(treat ~ pred, data,
      method = "multinom",
      subset = "subset_fail"
    ),
    regexp = "logical values"
  )

  expect_no_error(estimate_gps(treat ~ pred, data,
    method = "multinom",
    subset = "subset_pass"
  ))
})

## --testing methods and scenarios----------------------------------------------
test_that("estimate_gps: methods and scenarios", {
  data <- data.frame(
    treat = rep(c(1, 2, 3, 4, 5), 20),
    treat_fold = factor(rep(c(1, 2, 3, 4, 5), 20),
      ordered = TRUE
    ),
    treat_bin_log = rep(c(TRUE, FALSE), 50),
    treat_bin_char = rep(c("A", "B"), 50),
    pred = runif(100),
    subset_pass = rep(c(TRUE, FALSE), 50),
    subset_fail = rep(10, 100)
  )

  expect_no_error(estimate_gps(treat ~ pred, data, method = "multinom"))
  expect_no_error(estimate_gps(treat ~ pred, data, method = "vglm"))
  expect_no_error(estimate_gps(treat ~ pred, data, method = "brglm2"))
  expect_no_error(estimate_gps(treat ~ pred, data, method = "mblogit"))
  expect_error(estimate_gps(treat_fold ~ pred, data, method = "polr"),
    regexp = "ordered factor"
  )
})
