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

  ## testing
  expect_output(csregion(gps_matrix))
  expect_output(csregion(gps_matrix2))
  expect_error(csregion(data.frame()), regexp = "gps")

  # checking if borders work
  invisible(capture.output({
    s1 <- sum(attr(csregion(gps_matrix3, "include"), "filter_vector"))
    s2 <- sum(attr(csregion(gps_matrix3, "exclude"), "filter_vector"))
  }))

  expect_true(s1 > s2)
})
