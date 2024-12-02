#test the csregion() with default multinom() gps estimation method
test_that('match_gps checking arguments: csmatrix', {
  data <- data.frame(treat = rep(c(1, 2, 3, 4, 5), 20),
                     y = rep(c(TRUE, FALSE), 50),
                     pred = rnorm(100, 30, 8))

  gps_matrix <- estimate_gps(treat ~ pred, data, method = 'multinom')

  csmatrix <- csregion(gps_matrix)

  ## testing a clear run
  expect_no_error(invisible(match_gps(csmatrix)))

  ## testing class
  expect_error(match_gps(data), regexp = 'class')

  ## testing NULL
  expect_error(match_gps(NULL), regexp = 'missing')

  ## testing reference
  expect_no_error(match_gps(csmatrix, reference = '1'))
  expect_error(match_gps(csmatrix, reference = 'a'), regexp = 'unique')
  expect_error(match_gps(csmatrix, reference = FALSE), regexp = 'string')

  ## testing caliper
  expect_no_error(match_gps(csmatrix, caliper = 1))
  expect_error(match_gps(csmatrix, caliper = -1.1), regexp = 'positive')
  expect_error(match_gps(csmatrix, caliper = NULL), regexp = 'can not')
  expect_error(match_gps(csmatrix, caliper = 'a'),  regexp = 'numeric')
  expect_error(match_gps(csmatrix, caliper = c(1, 2, 3)),  regexp = 'length')
  expect_error(match_gps(csmatrix, caliper = rep('a', 10)),  regexp = 'numeric')
  expect_error(match_gps(csmatrix, caliper = rep(-1, 10)),  regexp = 'positive')
  expect_no_error(match_gps(csmatrix, caliper = rep(0.1, 10)))

  ## testing ratio
  expect_no_error(match_gps(csmatrix, ratio = 1))
  expect_no_error(match_gps(csmatrix, ratio = c(1:10)))
  expect_error(match_gps(csmatrix, ratio = c(1, 2)), regexp = 'atomic')
  expect_error(match_gps(csmatrix, ratio = rep('a', 10)), regexp = 'integer')
  expect_error(match_gps(csmatrix, ratio = 'a'), regexp = 'integer')
  expect_error(match_gps(csmatrix, ratio = 1.1), regexp = 'integer')

  ## testing replace
  expect_no_error(match_gps(csmatrix, replace = FALSE))
  expect_error(match_gps(csmatrix, replace = c(FALSE, TRUE)), regexp = 'length')
  expect_error(match_gps(csmatrix, replace = rep('a', 10)), regexp = 'logical')
  expect_no_error(match_gps(csmatrix, replace = rep(TRUE, 10)))
  expect_error(match_gps(csmatrix, replace = 'a'), regexp = 'flag')

  ## testing combos
  combos_fail1 <- c(1, 2, 3)
  combos_fail2 <- data.frame(a = c(1, 2, 3), b = c(1, 2, 3), c = c(1, 2, 3))
  combos_fail3 <- data.frame(a = c('a'), b = c('b'))
  combos_fail4 <- data.frame(a = c(1), b = c(1))
  combos_fail5 <- data.frame(a = c(1, 2), b = c(2, 1))
  combos_pass <- data.frame(a = c(1, 2, 3), b = c(2, 3, 5))

  expect_error(match_gps(csmatrix, combos = combos_fail1), regexp = 'data.frame')
  expect_error(match_gps(csmatrix, combos = combos_fail2), regexp = 'columns')
  expect_error(match_gps(csmatrix, combos = combos_fail3), regexp = 'unique')
  expect_error(match_gps(csmatrix, combos = combos_fail4), regexp = 'match')
  expect_error(match_gps(csmatrix, combos = combos_fail5), regexp = 'combination')
  expect_no_error(match_gps(csmatrix, combos = combos_pass))

  ## kmeans.args
  expect_no_error(match_gps(csmatrix, kmeans.args = list()))

  ## kmeans.cluster
  expect_error(match_gps(csmatrix, kmeans.cluster = NULL), regexp = 'NULL')
  expect_error(match_gps(csmatrix, kmeans.cluster = 'a'), regexp = 'integer')
  expect_error(match_gps(csmatrix, kmeans.cluster = rep('a', 10)), regexp = 'atomic')
  expect_error(match_gps(csmatrix, kmeans.cluster = c(1, 2, 3)), regexp = 'atomic')
  expect_error(match_gps(csmatrix, kmeans.cluster = rep(-1, 10)), regexp = 'greater')
  expect_error(match_gps(csmatrix, kmeans.cluster = rep(1, 10)), regexp = 'greater')
  expect_no_error(match_gps(csmatrix, kmeans.cluster = 4))
  expect_no_error(match_gps(csmatrix, kmeans.cluster = rep(4, 10)))
})
