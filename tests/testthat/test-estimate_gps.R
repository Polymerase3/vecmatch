## --testing formals: formula & data -------------------------------------------
test_that('Formals checking: formula', {
  data <- data.frame(y = runif(20),
                     group = rep(c(TRUE, FALSE), 10),
                     sex = rep(c('M', 'F', 'F', 'M'), 5))
  treat <- rep(c(0, 1), each = 10)
  treat_fail <- c(rep(c(0, 1), 9), c(NA, NA))
  pred <- runif(20)
  pred_fail <- runif(21)

  expect_error(estimate_gps(), regexp = 'missing')
  expect_error(estimate_gps(NULL), regexp = 'treatment')
  expect_error(estimate_gps(treat ~ pred_fail), regexp = 'samples')
  expect_error(estimate_gps(treat_fail ~ pred), regexp = 'NA')
  expect_warning(estimate_gps(data$y ~ data$group))
  expect_no_error(estimate_gps(treat ~ pred))
  expect_warning(estimate_gps(y ~ group * sex, data))
})

## --testing formals: method, ref, logicals------------------------------------------------------
test_that('Formals checking:  method, ref, logicals', {
  data <- data.frame(y = rep(c(TRUE, FALSE), 10),
                     pred = runif(20))

  #method
  expect_error(estimate_gps(y ~ pred, data, method = c()), regexp = 'string')
  expect_error(estimate_gps(y ~ pred, data, method = 'error'), regexp = 'method')
  expect_no_error(estimate_gps(y ~ pred, data, method = NULL))

  #ref
  expect_error(estimate_gps(y ~ pred, data, method = NULL, reference = TRUE),
               regexp = 'reference')
  expect_error(estimate_gps(y ~ pred, data, method = NULL, reference = "FAIL"),
               regexp = 'unique')
  expect_no_error(estimate_gps(y ~ pred, data, method = NULL, reference = NULL))

  #logicals
  expect_error(estimate_gps(y ~ pred, data, fit.object = 'fail'),
               regexp = 'logical')
  expect_error(estimate_gps(y ~ pred, data, verbose.output = 'fail'),
               regexp = 'logical')
  expect_no_error(estimate_gps(y ~ pred, data, verbose.output = TRUE))
})

## --testing formals: missing
test_that('Formals checking: missing and by', {
  data <- data.frame(treat = rep(c('A', 'B', 'C'), 7),
                       y = runif(21),
                     group = rep(c(TRUE, FALSE, TRUE), 7),
                     sex = c(rep(c('M', 'F', 'F', 'M'), 5), 'M'))

  #missing
  expect_error(estimate_gps(data$treat ~ data$y, missing = 'failure'),
               regexp = 'allowed')
  expect_error(estimate_gps(data$treat ~ data$y, missing = TRUE),
               regexp = 'string')
  expect_no_error(estimate_gps(data$treat ~ data$y,
                               missing = 'complete.cases'))
})

##--testing formals: ordinal.treat----------------------------------------------
test_that('Formals checking: ordinal.treat', {
  data <- data.frame(treat = rep(c(1, 2, 3, 4, 5), 20),
                     pred = runif(100))
  data$treat <- factor(data$treat, levels = c(1, 2, 3, 4, 5), ordered = TRUE)

  expect_error(estimate_gps(treat ~ pred, data, method = 'multinom',
                            ordinal.treat = list(1)), regexp = 'atomic')
  expect_error(estimate_gps(treat ~ pred, data, method = 'multinom',
                            ordinal.treat = c(1, 2, 3)), regexp = 'levels')
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'multinom',
                               ordinal.treat = c(1, 3, 2, 5, 4)))
})
