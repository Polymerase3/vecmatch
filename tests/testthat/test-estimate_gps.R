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

  ## reproducibility tests ----------------------------------------------------
  set.seed(1001)
  old_seed <- .Random.seed
  res1 <- estimate_gps(treat ~ pred)
  expect_identical(.Random.seed, old_seed)

  set.seed(2002)
  gps1 <- estimate_gps(treat ~ pred)
  set.seed(2002)
  gps2 <- estimate_gps(treat ~ pred)
  expect_identical(gps1, gps2)
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

  ## reproducibility tests ----------------------------------------------------
  set.seed(3003)
  old_seed <- .Random.seed
  res2 <- estimate_gps(y ~ pred, data, method = "multinom")
  expect_identical(.Random.seed, old_seed)

  set.seed(4004)
  g1 <- estimate_gps(y ~ pred, data, method = "multinom")
  set.seed(4004)
  g2 <- estimate_gps(y ~ pred, data, method = "multinom")
  expect_identical(g1, g2)
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

  ## reproducibility tests ----------------------------------------------------
  set.seed(5005)
  old_seed <- .Random.seed
  res_by <- estimate_gps(treat ~ y, data,
    method = "multinom",
    by = "sex"
  )
  expect_identical(.Random.seed, old_seed)

  set.seed(6006)
  g_by1 <- estimate_gps(treat ~ y, data,
    method = "multinom",
    by = "sex"
  )
  set.seed(6006)
  g_by2 <- estimate_gps(treat ~ y, data,
    method = "multinom",
    by = "sex"
  )
  expect_identical(g_by1, g_by2)
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

  ## reproducibility tests ----------------------------------------------------
  set.seed(7007)
  old_seed <- .Random.seed
  res_link <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    link = "generalized_logit"
  )
  expect_identical(.Random.seed, old_seed)

  set.seed(8008)
  g_link1 <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    link = "generalized_logit"
  )
  set.seed(8008)
  g_link2 <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    link = "generalized_logit"
  )
  expect_identical(g_link1, g_link2)
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

  ## reproducibility tests ----------------------------------------------------
  set.seed(9009)
  old_seed <- .Random.seed
  res_ord <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    ordinal_treat = c(1, 3, 2, 5, 4)
  )
  expect_identical(.Random.seed, old_seed)

  set.seed(10101)
  g_ord1 <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    ordinal_treat = c(1, 3, 2, 5, 4)
  )
  set.seed(10101)
  g_ord2 <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    ordinal_treat = c(1, 3, 2, 5, 4)
  )
  expect_identical(g_ord1, g_ord2)
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

  ## reproducibility tests ----------------------------------------------------
  set.seed(11111)
  old_seed <- .Random.seed
  res_sub <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    subset = "subset_pass"
  )
  expect_identical(.Random.seed, old_seed)

  set.seed(12121)
  g_sub1 <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    subset = "subset_pass"
  )
  set.seed(12121)
  g_sub2 <- estimate_gps(treat ~ pred, data,
    method = "multinom",
    subset = "subset_pass"
  )
  expect_identical(g_sub1, g_sub2)
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

  ## reproducibility tests ----------------------------------------------------
  set.seed(13131)
  old_seed <- .Random.seed
  res_meth <- estimate_gps(treat ~ pred, data, method = "multinom")
  expect_identical(.Random.seed, old_seed)

  set.seed(14141)
  g_meth1 <- estimate_gps(treat ~ pred, data, method = "brglm2")
  set.seed(14141)
  g_meth2 <- estimate_gps(treat ~ pred, data, method = "brglm2")
  expect_identical(g_meth1, g_meth2)
})

test_that("gps methods: setup helper object", {
  # small helper to build a minimal gps object
  make_test_gps <- function() {
    df <- data.frame(
      treatment = factor(c("A", "A", "B", "B")),
      A         = c(0.6, 0.4, 0.2, 0.1),
      B         = c(0.4, 0.6, 0.8, 0.9)
    )

    original_data <- data.frame(
      status = c("A", "A", "B", "B"),
      x      = 1:4
    )

    structure(
      df,
      original_data = original_data,
      function_call = quote(estimate_gps(
        formula = status ~ x,
        data = original_data
      )),
      class = c("gps", "data.frame")
    )
  }

  gps_obj <- make_test_gps()

  expect_s3_class(gps_obj, "gps")
  expect_true("treatment" %in% names(gps_obj))
  expect_equal(nrow(gps_obj), 4L)
  expect_equal(ncol(gps_obj), 3L)

  ## reproducibility tests ----------------------------------------------------
  set.seed(15151)
  old_seed <- .Random.seed
  tmp_obj <- make_test_gps()
  expect_identical(.Random.seed, old_seed)

  set.seed(16161)
  g_make1 <- make_test_gps()
  set.seed(16161)
  g_make2 <- make_test_gps()
  expect_identical(g_make1, g_make2)
})

test_that(".print_gps_core prints header and returns invisibly", {
  make_test_gps <- function() {
    df <- data.frame(
      treatment = factor(c("A", "A", "B", "B")),
      A         = c(0.6, 0.4, 0.2, 0.1),
      B         = c(0.4, 0.6, 0.8, 0.9)
    )
    structure(df, class = c("gps", "data.frame"))
  }

  gps_obj <- make_test_gps()

  utils::capture.output({
    out <- utils::capture.output(
      {
        res <- .print_gps_core(gps_obj)
      },
      type = "message"
    )
  })

  # check that the header line is present
  expect_true(any(grepl("number of units: 4", tolower(out), fixed = TRUE)))

  # helper returns the object invisibly
  expect_identical(res, gps_obj)

  # class and structure are preserved
  expect_s3_class(res, "gps")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 4L)
  expect_equal(colnames(res), c("treatment", "A", "B"))

  ## reproducibility tests ----------------------------------------------------
  set.seed(17171)
  old_seed <- .Random.seed
  utils::capture.output(.print_gps_core(gps_obj), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(18181)
  out1 <- utils::capture.output(.print_gps_core(gps_obj), type = "message")
  set.seed(18181)
  out2 <- utils::capture.output(.print_gps_core(gps_obj), type = "message")
  expect_identical(out1, out2)
})

test_that("print.gps uses .print_gps_core and prints header", {
  make_test_gps <- function() {
    df <- data.frame(
      treatment = factor(c("A", "A", "B", "B")),
      A         = c(0.6, 0.4, 0.2, 0.1),
      B         = c(0.4, 0.6, 0.8, 0.9)
    )
    structure(df, class = c("gps", "data.frame"))
  }

  gps_obj <- make_test_gps()

  utils::capture.output({
    out <- utils::capture.output(
      {
        res <- print(gps_obj)
      },
      type = "message"
    )
  })

  # header from print.gps
  expect_true(any(grepl("gps object \\(generalized propensity scores\\)", out)))

  expect_identical(res, gps_obj)

  ## reproducibility tests ----------------------------------------------------
  set.seed(19191)
  old_seed <- .Random.seed
  utils::capture.output(print(gps_obj), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(20202)
  out3 <- utils::capture.output(print(gps_obj), type = "message")
  set.seed(20202)
  out4 <- utils::capture.output(print(gps_obj), type = "message")
  expect_identical(out3, out4)
})

test_that("str.gps prints structured header and underlying structure", {
  make_test_gps <- function() {
    df <- data.frame(
      treatment = factor(c("A", "A", "B", "B")),
      A         = c(0.6, 0.4, 0.2, 0.1),
      B         = c(0.4, 0.6, 0.8, 0.9)
    )
    od <- data.frame(status = c("A", "A", "B", "B"), x = 1:4)
    structure(
      df,
      original_data = od,
      function_call = quote(estimate_gps(formula = status ~ x, data = od)),
      class         = c("gps", "data.frame")
    )
  }

  gps_obj <- make_test_gps()

  out <- testthat::capture_output({
    res <- str(gps_obj)
  })

  # header lines
  expect_true(grepl("gps object: generalized propensity scores", out, fixed = TRUE))
  expect_true(grepl("Dimensions: 4 rows x 3 columns", out, fixed = TRUE))
  expect_true(grepl("Treatment column: treatment", out, fixed = TRUE))
  expect_true(grepl("GPS columns: A, B", out, fixed = TRUE))

  # underlying structure mention
  expect_true(grepl("Underlying data.frame structure:", out, fixed = TRUE))

  expect_identical(res, gps_obj)

  ## reproducibility tests ----------------------------------------------------
  set.seed(21212)
  old_seed <- .Random.seed
  testthat::capture_output(str(gps_obj))
  expect_identical(.Random.seed, old_seed)

  set.seed(22222)
  out_s1 <- testthat::capture_output(str(gps_obj))
  set.seed(22222)
  out_s2 <- testthat::capture_output(str(gps_obj))
  expect_identical(out_s1, out_s2)
})

test_that("summary.gps returns structured summary.gps object", {
  make_test_gps <- function() {
    df <- data.frame(
      treatment = factor(c("A", "A", "B", "B")),
      A         = c(0.6, 0.4, 0.2, 0.1),
      B         = c(0.4, 0.6, 0.8, 0.9)
    )
    od <- data.frame(status = c("A", "A", "B", "B"), x = 1:4)
    structure(
      df,
      original_data = od,
      function_call = quote(estimate_gps(formula = status ~ x, data = od)),
      class         = c("gps", "data.frame")
    )
  }

  gps_obj <- make_test_gps()

  s <- summary(gps_obj)

  expect_s3_class(s, "summary.gps")
  expect_equal(s$n, 4L)
  expect_equal(s$treatments, c("A", "B"))
  expect_true(is.list(s$gps_by_treatment))
  expect_setequal(names(s$gps_by_treatment), c("A", "B"))

  # each element of gps_by_treatment should be a matrix with rows = summary stats
  expect_true(all(vapply(s$gps_by_treatment, is.matrix, logical(1L))))

  # original_data carried through
  expect_true(is.data.frame(s$original_data))
  expect_equal(nrow(s$original_data), 4L)

  ## reproducibility tests ----------------------------------------------------
  set.seed(23232)
  old_seed <- .Random.seed
  s_tmp <- summary(gps_obj)
  expect_identical(.Random.seed, old_seed)

  set.seed(24242)
  s1 <- summary(gps_obj)
  set.seed(24242)
  s2 <- summary(gps_obj)
  expect_identical(s1, s2)
})

test_that("print.summary.gps prints without error and includes sections", {
  make_test_gps <- function() {
    df <- data.frame(
      treatment = factor(c("A", "A", "B", "B")),
      A         = c(0.6, 0.4, 0.2, 0.1),
      B         = c(0.4, 0.6, 0.8, 0.9)
    )
    od <- data.frame(status = c("A", "A", "B", "B"), x = 1:4)
    structure(
      df,
      original_data = od,
      function_call = quote(estimate_gps(formula = status ~ x, data = od)),
      class         = c("gps", "data.frame")
    )
  }

  gps_obj <- make_test_gps()
  s <- summary(gps_obj)

  out <- utils::capture.output(
    {
      res <- print(s)
    },
    type = "message"
  )

  out <- cli::ansi_strip(out)

  expect_true(any(grepl("Summary of gps object", out, fixed = TRUE)))
  expect_true(any(grepl("GPS matrix by treatment", out, fixed = TRUE)))
  expect_true(any(grepl("Original data summary", out, fixed = TRUE)))

  expect_identical(res, s)

  ## reproducibility tests ----------------------------------------------------
  set.seed(25252)
  old_seed <- .Random.seed
  utils::capture.output(print(s), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(26262)
  out_ps1 <- utils::capture.output(print(s), type = "message")
  set.seed(26262)
  out_ps2 <- utils::capture.output(print(s), type = "message")
  expect_identical(out_ps1, out_ps2)
})

test_that("as.data.frame.gps drops gps class and keeps data", {
  df <- data.frame(
    treatment = factor(c("A", "A", "B", "B")),
    A         = c(0.6, 0.4, 0.2, 0.1),
    B         = c(0.4, 0.6, 0.8, 0.9)
  )

  gps_obj <- structure(df, class = c("gps", "data.frame"))

  df2 <- as.data.frame(gps_obj)

  expect_s3_class(df2, "data.frame")
  expect_false(inherits(df2, "gps"))
  expect_identical(df2, df)

  ## reproducibility tests ----------------------------------------------------
  set.seed(27272)
  old_seed <- .Random.seed
  tmp_df <- as.data.frame(gps_obj)
  expect_identical(.Random.seed, old_seed)

  set.seed(28282)
  df_a <- as.data.frame(gps_obj)
  set.seed(28282)
  df_b <- as.data.frame(gps_obj)
  expect_identical(df_a, df_b)
})

test_that("plot.gps runs without error and returns invisibly", {
  df <- data.frame(
    treatment = factor(rep(c("A", "B"), each = 5)),
    A         = seq(0.1, 1.0, length.out = 10),
    B         = seq(1.0, 0.1, length.out = 10)
  )
  gps_obj <- structure(df, class = c("gps", "data.frame"))

  # should not error and should return the input object invisibly
  ret <- NULL
  expect_silent({
    ret <- plot(gps_obj)
  })
  expect_identical(ret, gps_obj)

  ## reproducibility tests ----------------------------------------------------
  set.seed(29292)
  old_seed <- .Random.seed
  utils::capture.output(plot(gps_obj), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(30303)
  r1 <- plot(gps_obj)
  set.seed(30303)
  r2 <- plot(gps_obj)
  expect_identical(r1, r2)
})

test_that("plot.gps errors when gps_cols are invalid", {
  df <- data.frame(
    treatment = factor(rep(c("A", "B"), each = 3)),
    A         = seq(0.1, 0.6, length.out = 6),
    B         = seq(0.6, 0.1, length.out = 6)
  )
  gps_obj <- structure(df, class = c("gps", "data.frame"))

  expect_error(
    plot(gps_obj, gps_cols = "C"),
    "None of the requested `gps_cols` are valid",
    fixed = TRUE
  )

  ## reproducibility tests ----------------------------------------------------
  set.seed(31313)
  old_seed <- .Random.seed
  expect_error(
    plot(gps_obj, gps_cols = "C"),
    "None of the requested `gps_cols` are valid",
    fixed = TRUE
  )
  expect_identical(.Random.seed, old_seed)
})
