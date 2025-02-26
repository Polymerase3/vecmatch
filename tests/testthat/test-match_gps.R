# test the csregion() with default multinom() gps estimation method
test_that("match_gps checking arguments: csmatrix", {
  withr::with_seed(1643741, {
    data <- data.frame(
      treat = rep(c(1, 2, 3, 4, 5), 60),
      y = rep(c(TRUE, FALSE), 150),
      pred = rnorm(300, 30, 8)
    )
  })

  # estimate the gps
  gps_matrix <- estimate_gps(treat ~ pred, data, method = "multinom")

  # drop observations outside the csr
  invisible(capture.output(
    {
      csmatrix <- csregion(gps_matrix)
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
    expect_no_error(match_gps(csmatrix, reference = "1"))
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

    ## testing combos
    combos_fail1 <- c(1, 2, 3)
    combos_fail2 <- data.frame(a = c(1, 2, 3), b = c(1, 2, 3), c = c(1, 2, 3))
    combos_fail3 <- data.frame(a = c("a"), b = c("b"))
    combos_fail4 <- data.frame(a = c(1), b = c(1))
    combos_fail5 <- data.frame(a = c(1, 2), b = c(2, 1))
    combos_pass <- data.frame(a = c(1, 2, 3), b = c(2, 3, 5))

    expect_error(match_gps(csmatrix, combos = combos_fail1),
      regexp = "data.frame"
    )
    expect_error(match_gps(csmatrix, combos = combos_fail2),
      regexp = "columns"
    )
    expect_error(match_gps(csmatrix, combos = combos_fail3),
      regexp = "unique"
    )
    expect_error(match_gps(csmatrix, combos = combos_fail4),
      regexp = "match"
    )
    expect_error(match_gps(csmatrix, combos = combos_fail5),
      regexp = "combination"
    )

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
    expect_no_error(match_gps(csmatrix, reference = "1", method = "fullopt"))
  })
})
