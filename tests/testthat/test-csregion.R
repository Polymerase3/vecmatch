# test the csregion() with default multinom() gps estimation method
test_that("csregion() with multinom data", {
  withr::with_seed(6134423, {
    data <- data.frame(
      treat = rep(c(1, 2, 3, 4, 5), 120),
      y = rep(c(TRUE, FALSE), 300),
      pred = runif(600),
      pred2 = rnorm(600)
    )
  })

  gps_matrix <- estimate_gps(treat ~ pred, data, method = "multinom")
  gps_matrix2 <- estimate_gps(y ~ pred, data, method = "multinom")
  gps_matrix3 <- estimate_gps(treat ~ pred2, data)

  ## testing -------------------------------------------------------------------
  expect_error(csregion(data.frame()), regexp = "gps")

  # checking if borders work
  s1 <- sum(attr(csregion(gps_matrix3, "include"), "filter_vector"))
  s2 <- sum(attr(csregion(gps_matrix3, "exclude"), "filter_vector"))

  expect_true(s1 > s2)

  ## reproducibility tests ----------------------------------------------------
  # csregion() should not change the global .Random.seed
  set.seed(12345)
  old_seed <- .Random.seed
  tmp_csr <- csregion(gps_matrix3)
  expect_identical(.Random.seed, old_seed)

  # running csregion() twice with the same seed should give identical results
  set.seed(4242)
  csr1 <- csregion(gps_matrix3)
  set.seed(4242)
  csr2 <- csregion(gps_matrix3)
  expect_identical(csr1, csr2)
})

test_that("csr methods: summary/print/str/plot/as.data.frame", {
  # build a small csr object via the usual pipeline
  withr::with_seed(6134423, {
    data <- data.frame(
      treat = rep(c(1, 2, 3, 4, 5), 120),
      y     = rep(c(TRUE, FALSE), 300),
      pred  = runif(600),
      pred2 = rnorm(600)
    )
  })

  gps_matrix <- estimate_gps(treat ~ pred, data, method = "multinom")
  csr_obj <- csregion(gps_matrix)

  expect_s3_class(csr_obj, "csr")
  expect_s3_class(csr_obj, "gps")
  expect_s3_class(csr_obj, "data.frame")

  ## summary.csr ---------------------------------------------------------------
  s <- summary(csr_obj)

  expect_s3_class(s, "summary.csr")
  expect_true(is.numeric(s$n_total))
  expect_true(is.numeric(s$n_kept))
  expect_true(s$n_total >= s$n_kept)
  expect_true(is.data.frame(s$bounds_table))
  expect_true(is.data.frame(s$groups_table))
  expect_true(all(c("Treatment", "n_before", "n_after", "Excluded") %in%
    names(s$groups_table)))

  ## print.summary.csr ---------------------------------------------------------
  out_sum <- utils::capture.output(
    {
      res_sum <- print(s)
    },
    type = "message"
  )

  out_sum_clean <- cli::ansi_strip(out_sum)
  expect_true(any(grepl("Rectangular CSR Borders Evaluation", out_sum_clean,
    fixed = TRUE
  )))

  expect_identical(res_sum, s)

  ## print.csr -----------------------------------------------------------------
  out_pc <- utils::capture.output(
    {
      res_pc <- print(csr_obj)
    },
    type = "message"
  )

  out_pc_clean <- cli::ansi_strip(out_pc)
  expect_true(any(grepl(
    "csr object \\(gps filtered to common support region\\)",
    out_pc_clean
  )))

  expect_identical(res_pc, csr_obj)

  ## str.csr -------------------------------------------------------------------
  out_sc <- utils::capture.output({
    res_sc <- str(csr_obj)
  })

  expect_true(any(grepl("csr object: gps filtered to common support region",
    out_sc,
    fixed = TRUE
  )))
  expect_identical(res_sc, csr_obj)

  ## as.data.frame.csr ---------------------------------------------------------
  df_csr <- as.data.frame(csr_obj)

  expect_s3_class(df_csr, "data.frame")
  expect_false(inherits(df_csr, "csr"))
  expect_equal(nrow(df_csr), nrow(csr_obj))
  expect_equal(colnames(df_csr), colnames(csr_obj))

  ## plot.csr ------------------------------------------------------------------
  res_plot <- NULL
  expect_silent({
    res_plot <- plot(csr_obj)
  })
  expect_identical(res_plot, csr_obj)

  ## reproducibility tests ----------------------------------------------------
  # summary.csr should not change global seed and should be reproducible
  set.seed(111)
  old_seed <- .Random.seed
  s_tmp <- summary(csr_obj)
  expect_identical(.Random.seed, old_seed)

  set.seed(222)
  s1 <- summary(csr_obj)
  set.seed(222)
  s2 <- summary(csr_obj)
  expect_identical(s1, s2)

  # print.csr should not change global seed
  set.seed(333)
  old_seed <- .Random.seed
  utils::capture.output(print(csr_obj), type = "message")
  expect_identical(.Random.seed, old_seed)

  # str.csr should not change global seed
  set.seed(444)
  old_seed <- .Random.seed
  utils::capture.output(str(csr_obj))
  expect_identical(.Random.seed, old_seed)

  # as.data.frame.csr should be reproducible for a given seed
  set.seed(555)
  df1 <- as.data.frame(csr_obj)
  set.seed(555)
  df2 <- as.data.frame(csr_obj)
  expect_identical(df1, df2)

  # plot.csr should not change global seed and should be reproducible
  set.seed(666)
  old_seed <- .Random.seed
  utils::capture.output(plot(csr_obj), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(777)
  p1 <- plot(csr_obj)
  set.seed(777)
  p2 <- plot(csr_obj)
  expect_identical(p1, p2)
})

test_that("plot.csr handles gps_cols argument and layout branches", {
  withr::with_seed(6134423, {
    data <- data.frame(
      treat = rep(c(1, 2, 3, 4, 5), 120),
      y     = rep(c(TRUE, FALSE), 300),
      pred  = runif(600)
    )
  })

  gps_matrix <- estimate_gps(treat ~ pred, data, method = "multinom")
  csr_obj <- csregion(gps_matrix)

  all_gps_cols <- setdiff(names(csr_obj), "treatment")

  # sanity: we need at least 5 gps columns to hit the last branch
  expect_gte(length(all_gps_cols), 5L)

  # 1) default: gps_cols = NULL -> uses all_gps_cols
  expect_silent({
    res_default <- plot(csr_obj)
  })
  expect_identical(res_default, csr_obj)

  # 2) one gps column -> branch n_gps == 1
  expect_silent({
    res_1 <- plot(csr_obj, gps_cols = all_gps_cols[1])
  })
  expect_identical(res_1, csr_obj)

  # 3) two gps columns -> branch n_gps == 2
  expect_silent({
    res_2 <- plot(csr_obj, gps_cols = all_gps_cols[1:2])
  })
  expect_identical(res_2, csr_obj)

  # 4) three gps columns -> branch n_gps <= 4
  expect_silent({
    res_3 <- plot(csr_obj, gps_cols = all_gps_cols[1:3])
  })
  expect_identical(res_3, csr_obj)

  # 5) five gps columns -> branch n_gps > 4 (uses ceiling(n_gps / 2))
  expect_silent({
    res_5 <- plot(csr_obj, gps_cols = all_gps_cols[1:5])
  })
  expect_identical(res_5, csr_obj)

  # 6) invalid gps_cols -> error as expected
  expect_error(
    plot(csr_obj, gps_cols = "not_a_valid_column"),
    regexp = "None of the requested `gps_cols` are valid\\."
  )

  ## reproducibility tests ----------------------------------------------------
  # plot.csr with specific gps_cols should not change global seed
  set.seed(999)
  old_seed <- .Random.seed
  utils::capture.output(
    plot(csr_obj, gps_cols = all_gps_cols[1:5]),
    type = "message"
  )
  expect_identical(.Random.seed, old_seed)

  # running plot.csr twice with the same seed should return the same object
  set.seed(1010)
  p_a <- plot(csr_obj, gps_cols = all_gps_cols[1:4])
  set.seed(1010)
  p_b <- plot(csr_obj, gps_cols = all_gps_cols[1:4])
  expect_identical(p_a, p_b)
})
