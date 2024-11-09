#test the csregion() with default multinom() gps estimation method
test_that('csregion() with multinom data', {
  data <- data.frame(treat = rep(c(1, 2, 3, 4, 5), 20),
                     y = rep(c(TRUE, FALSE), 10),
                     pred = runif(20))

  gps_matrix <- estimate_gps(treat ~ pred, data, method = 'multinom')

  ## testing
  expect_no_error(invisible(csregion(gps_matrix)))
  expect_error(csregion(data.frame()), regexp = 'gps')
})

## fix the CSR for two treats
## add tratment column to all fitted_obj
