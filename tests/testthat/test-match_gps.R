# test the csregion() with default multinom() gps estimation method
test_that("match_gps checking arguments: csmatrix", {
  withr::with_seed(1643741, {
    data <- data.frame(
      treat = rep(c("A", "B", "C", "D", "E"), 60),
      y = rep(c(TRUE, FALSE), 150),
      pred = rnorm(300, 30, 8)
    )
  })

  # estimate the gps
  gps_matrix <- estimate_gps(treat ~ pred, data, method = "multinom")
  gps_matrix_2t <- estimate_gps(
    treat ~ pred,
    data[data$treat %in% c("A", "B"), ],
    method = "multinom"
  )

  # drop observations outside the csr
  invisible(capture.output(
    {
      csmatrix <- csregion(gps_matrix)
      csmatrix_2t <- csregion(gps_matrix_2t)
    },
    file = NULL
  ))

  ## testing a clear run
  withr::with_options(list(warn = -1), {
    expect_no_error(match_gps(csmatrix))

    ## testing class
    expect_error(match_gps(data), regexp = "class")

    ## testing NULL
    expect_error(match_gps(NULL), regexp = "missing")

    ## testing reference
    expect_no_error(match_gps(csmatrix, reference = "A"))
    expect_error(match_gps(csmatrix, reference = "a"), regexp = "unique")
    expect_error(match_gps(csmatrix, reference = FALSE), regexp = "string")

    ## testing caliper
    expect_no_error(match_gps(csmatrix, caliper = 1))
    expect_error(match_gps(csmatrix, caliper = -1.1), regexp = "positive")
    expect_error(match_gps(csmatrix, caliper = "a"), regexp = "numeric")
    expect_error(match_gps(csmatrix, caliper = c(1, 2, 3)), regexp = "length")

    ## testing ratio
    expect_no_error(match_gps(csmatrix, ratio = 1))
    expect_error(match_gps(csmatrix, ratio = c(1:4)), regexp = "matches")
    expect_error(match_gps(csmatrix, ratio = c(1, 2)), regexp = "atomic")
    expect_error(match_gps(csmatrix, ratio = rep("a", 10)), regexp = "integer")

    ## testing replace
    expect_no_error(match_gps(csmatrix, replace = FALSE))
    expect_error(
      match_gps(csmatrix, replace = c(FALSE, TRUE)),
      regexp = "length"
    )
    expect_error(
      match_gps(csmatrix, replace = rep("a", 10)),
      regexp = "logical"
    )
    expect_no_error(match_gps(csmatrix, replace = rep(TRUE, 4)))

    ## kmeans.args
    expect_no_error(match_gps(csmatrix, kmeans.args = list()))

    ## kmeans_cluster
    expect_error(match_gps(csmatrix, kmeans_cluster = NULL), regexp = "NULL")
    expect_error(match_gps(csmatrix, kmeans_cluster = "a"), regexp = "integer")
    expect_error(
      match_gps(csmatrix, kmeans_cluster = rep("a", 10)),
      regexp = "atomic"
    )
    expect_error(match_gps(csmatrix, kmeans_cluster = -1), regexp = "greater")
    expect_error(
      match_gps(csmatrix, kmeans_cluster = rep(1, 10)),
      regexp = "equal"
    )
    expect_no_error(match_gps(csmatrix, kmeans_cluster = 4))
    expect_no_error(match_gps(csmatrix, kmeans_cluster = rep(4, 4)))

    ## matching methods
    expect_no_error(match_gps(csmatrix, reference = "A", method = "fullopt"))

    ## test the two treatments cases for "nnm" and "fullopt"
    expect_no_error(match_gps(csmatrix_2t, reference = "A", method = "nnm"))
    expect_no_error(match_gps(csmatrix_2t, reference = "A", method = "fullopt"))
  })

  ## reproducibility tests ----------------------------------------------------
  # match_gps should not change the global RNG state
  set.seed(424242)
  old_seed <- .Random.seed
  tmp_match <- match_gps(
    csmatrix,
    reference      = "A",
    caliper        = 1,
    kmeans_cluster = 4
  )
  expect_identical(.Random.seed, old_seed)

  # running match_gps twice with the same seed should give identical results
  set.seed(434343)
  match1 <- match_gps(
    csmatrix,
    reference      = "A",
    caliper        = 1,
    kmeans_cluster = 4
  )
  set.seed(434343)
  match2 <- match_gps(
    csmatrix,
    reference      = "A",
    caliper        = 1,
    kmeans_cluster = 4
  )
  expect_identical(match1, match2)
})

test_that("matched methods: str/print/summary/plot/as.data.frame", {
  # build a small matched object via the usual pipeline
  withr::with_seed(1643741, {
    data <- data.frame(
      treat = rep(c("A", "B", "C", "D", "E"), 60),
      y     = rep(c(TRUE, FALSE), 150),
      pred  = rnorm(300, 30, 8)
    )
  })

  gps_matrix <- estimate_gps(treat ~ pred, data, method = "multinom")

  invisible(capture.output(
    {
      csmatrix <- csregion(gps_matrix)
    },
    file = NULL
  ))

  matched_obj <- match_gps(
    csmatrix,
    reference      = "A",
    caliper        = 1,
    kmeans_cluster = 4
  )

  expect_s3_class(matched_obj, "matched")
  expect_s3_class(matched_obj, "data.frame")

  ## str.matched ---------------------------------------------------------------
  out_str <- utils::capture.output({
    res_str <- str(matched_obj)
  })

  expect_true(any(grepl("Matched data.frame:", out_str, fixed = TRUE)))
  expect_identical(res_str, matched_obj)

  ## .print_matched_core -------------------------------------------------------
  out_core <- utils::capture.output(
    {
      res_core <- .print_matched_core(matched_obj)
    },
    type = "message"
  )

  # strip ansi if present and check header line
  out_core_clean <- cli::ansi_strip(out_core)
  expect_true(any(grepl("Number of units after matching", out_core_clean,
    fixed = TRUE
  )))

  expect_identical(res_core, matched_obj)

  ## print.matched -------------------------------------------------------------
  out_print <- utils::capture.output(
    {
      res_print <- print(matched_obj)
    },
    type = "message"
  )

  out_print_clean <- cli::ansi_strip(out_print)
  expect_true(any(grepl(
    "matched object \\(GPS-based matched dataset\\)",
    out_print_clean
  )))

  expect_identical(res_print, matched_obj)

  ## summary.matched -----------------------------------------------------------
  s <- summary(matched_obj)

  expect_s3_class(s, "summary.matched")
  expect_true(is.numeric(s$n_total))
  expect_true(is.numeric(s$n_matched))
  expect_true(s$n_total >= s$n_matched)
  expect_true(is.data.frame(s$per_treatment))
  expect_true(all(c("Treatment", "n_before", "n_after", "Retained_percent") %in%
    names(s$per_treatment)))
  expect_true(s$treatment_var %in% names(matched_obj))

  ## print.summary.matched -----------------------------------------------------
  out_ps <- utils::capture.output(
    {
      res_ps <- print(s)
    },
    type = "message"
  )

  out_ps_clean <- cli::ansi_strip(out_ps)
  expect_true(any(grepl("Summary of matched object", out_ps_clean,
    fixed = TRUE
  )))

  expect_identical(res_ps, s)

  ## plot.matched --------------------------------------------------------------
  res_plot <- NULL
  expect_silent({
    res_plot <- plot(matched_obj)
  })
  expect_identical(res_plot, matched_obj)

  ## as.data.frame.matched -----------------------------------------------------
  df_m <- as.data.frame(matched_obj)

  expect_s3_class(df_m, "data.frame")
  expect_false(inherits(df_m, "matched"))
  expect_equal(nrow(df_m), nrow(matched_obj))
  expect_equal(colnames(df_m), colnames(matched_obj))

  ## reproducibility tests -----------------------------------------------------
  # str.matched should not touch RNG and should be reproducible
  set.seed(515151)
  old_seed <- .Random.seed
  out_str2 <- utils::capture.output(str(matched_obj))
  expect_identical(.Random.seed, old_seed)

  set.seed(525252)
  out_str_a <- utils::capture.output(str(matched_obj))
  set.seed(525252)
  out_str_b <- utils::capture.output(str(matched_obj))
  expect_identical(out_str_a, out_str_b)

  # .print_matched_core should not touch RNG and be reproducible
  set.seed(535353)
  old_seed <- .Random.seed
  out_core2 <- utils::capture.output(.print_matched_core(matched_obj),
    type = "message"
  )
  expect_identical(.Random.seed, old_seed)

  set.seed(545454)
  out_core_a <- utils::capture.output(.print_matched_core(matched_obj),
    type = "message"
  )
  set.seed(545454)
  out_core_b <- utils::capture.output(.print_matched_core(matched_obj),
    type = "message"
  )
  expect_identical(out_core_a, out_core_b)

  # print.matched should not touch RNG and be reproducible
  set.seed(555555)
  old_seed <- .Random.seed
  out_print2 <- utils::capture.output(print(matched_obj), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(565656)
  out_print_a <- utils::capture.output(print(matched_obj), type = "message")
  set.seed(565656)
  out_print_b <- utils::capture.output(print(matched_obj), type = "message")
  expect_identical(out_print_a, out_print_b)

  # summary.matched should not touch RNG and be reproducible
  set.seed(575757)
  old_seed <- .Random.seed
  s2 <- summary(matched_obj)
  expect_identical(.Random.seed, old_seed)

  set.seed(585858)
  s_a <- summary(matched_obj)
  set.seed(585858)
  s_b <- summary(matched_obj)
  expect_identical(s_a, s_b)

  # print.summary.matched should not touch RNG and be reproducible
  set.seed(595959)
  old_seed <- .Random.seed
  out_ps2 <- utils::capture.output(print(s), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(606060)
  out_ps_a <- utils::capture.output(print(s), type = "message")
  set.seed(606060)
  out_ps_b <- utils::capture.output(print(s), type = "message")
  expect_identical(out_ps_a, out_ps_b)

  # plot.matched should not touch RNG and be reproducible w.r.t. returned object
  set.seed(616161)
  old_seed <- .Random.seed
  res_plot2 <- plot(matched_obj)
  expect_identical(.Random.seed, old_seed)

  set.seed(626262)
  p1 <- plot(matched_obj)
  set.seed(626262)
  p2 <- plot(matched_obj)
  expect_identical(p1, p2)

  # plot.matched() does not permanently change par settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  plot(matched_obj)

  new_par <- par(no.readonly = TRUE)
  expect_identical(old_par, new_par)

  # as.data.frame.matched should not touch RNG and be reproducible
  set.seed(636363)
  old_seed <- .Random.seed
  df_m2 <- as.data.frame(matched_obj)
  expect_identical(.Random.seed, old_seed)

  set.seed(646464)
  df_a <- as.data.frame(matched_obj)
  set.seed(646464)
  df_b <- as.data.frame(matched_obj)
  expect_identical(df_a, df_b)
})
