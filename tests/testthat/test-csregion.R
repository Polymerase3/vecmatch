# test the csregion() with default multinom() gps estimation method
test_that("csregion() with multinom data", {
  withr::with_seed(6134423, {
    data <- data.frame(
      treat = rep(c(1, 2, 3, 4, 5), 120),
      y = rep(c(TRUE, FALSE), 300),
      pred = runif(600)
    )
  })

  gps_matrix <- estimate_gps(treat ~ pred, data, method = "multinom")
  gps_matrix2 <- estimate_gps(y ~ pred, data, method = "multinom")

  ## testing
  expect_output(csregion(gps_matrix))
  expect_output(csregion(gps_matrix2))
  expect_error(csregion(data.frame()), regexp = "gps")
})
