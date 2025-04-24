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
  gps_matrix_2t <- estimate_gps(treat ~ pred, data[data$treat %in% c("A", "B"), ],
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
    expect_error(match_gps(csmatrix, replace = c(FALSE, TRUE)),
      regexp = "length"
    )
    expect_error(match_gps(csmatrix, replace = rep("a", 10)),
      regexp = "logical"
    )
    expect_no_error(match_gps(csmatrix, replace = rep(TRUE, 4)))

    ## kmeans.args
    expect_no_error(match_gps(csmatrix, kmeans.args = list()))

    ## kmeans_cluster
    expect_error(match_gps(csmatrix, kmeans_cluster = NULL), regexp = "NULL")
    expect_error(match_gps(csmatrix, kmeans_cluster = "a"), regexp = "integer")
    expect_error(match_gps(csmatrix, kmeans_cluster = rep("a", 10)),
      regexp = "atomic"
    )
    expect_error(match_gps(csmatrix, kmeans_cluster = -1), regexp = "greater")
    expect_error(match_gps(csmatrix, kmeans_cluster = rep(1, 10)),
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
})
