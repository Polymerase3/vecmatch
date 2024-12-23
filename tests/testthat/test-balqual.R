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
})
