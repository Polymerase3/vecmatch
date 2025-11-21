test_that("balqual: argument checks and basic run", {
  # estimate the gps
  gps_matrix <- estimate_gps(status ~ age,
    cancer,
    method = "multinom",
    refernce = "control"
  )

  # drop observations outside the csr
  invisible(capture.output(
    {
      csmatrix <- csregion(gps_matrix)
    },
    file = NULL
  ))

  ## testing a clear run
  withr::with_options(list(warn = -1), {
    # matching the csmatrix
    matched_cancer <- match_gps(csmatrix,
      reference = "control",
      caliper = 1,
      kmeans_cluster = 5
    )
  })

  # basic test run
  expect_no_error(
    invisible(
      capture.output(
        balqual(matched_cancer, status ~ age),
        file = NULL
      )
    )
  )

  # test mean
  expect_no_error(
    invisible(
      capture.output(
        balqual(matched_cancer, status ~ age, statistic = "mean"),
        file = NULL
      )
    )
  )

  # test max
  expect_no_error(
    invisible(
      capture.output(
        balqual(matched_cancer, status ~ age * sex, statistic = "max"),
        file = NULL
      )
    )
  )

  # break cutoffs
  expect_error(
    invisible(
      capture.output(
        balqual(matched_cancer, status ~ age, cutoffs = 1),
        file = NULL
      )
    ),
    regexp = "length"
  )

  # run balqual once to obtain a quality object
  quality_obj <- balqual(
    matched_cancer,
    status ~ age,
    statistic = "max"
  )

  # capture printed output from str()
  out <- utils::capture.output({
    res <- str(quality_obj)
  })

  # check header string is present somewhere in the output
  expect_true(any(grepl("quality object: matching diagnostics", out, fixed = TRUE)))

  # str() should return the object invisibly
  expect_identical(res, quality_obj)

  # basic structural sanity checks
  expect_s3_class(quality_obj, "quality")
  expect_true("type" %in% names(quality_obj))
  expect_true("statistic" %in% names(quality_obj))

  # attributes used by str.quality should exist (at least some of them)
  expect_true(!is.null(attr(quality_obj, "original_data_before")) ||
    !is.null(attr(quality_obj, "original_data_after")))

  ## reproducibility tests ----------------------------------------------------
  # 1) balqual() should not change the global .Random.seed
  set.seed(12345)
  old_seed <- .Random.seed

  tmp_quality <- balqual(
    matched_cancer,
    status ~ age,
    statistic = "max"
  )

  expect_identical(.Random.seed, old_seed)

  # 2) running balqual() twice with the same seed should give identical results
  set.seed(4242)
  q1 <- balqual(
    matched_cancer,
    status ~ age,
    statistic = "max"
  )

  set.seed(4242)
  q2 <- balqual(
    matched_cancer,
    status ~ age,
    statistic = "max"
  )

  expect_identical(q1, q2)
})
