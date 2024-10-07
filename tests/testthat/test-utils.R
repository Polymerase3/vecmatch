#--testing get_formula_vars()---------------------------------------------------
##data
test_that('Testing data argument', {
  data <- data.frame()
  datax <- data.frame(treat = runif(20),
                      pred = rep(c(TRUE, FALSE), 10))
  expect_no_error(.get_formula_vars(datax$treat ~ datax$pred))
  expect_no_error(.get_formula_vars(treat ~ pred, data = datax))
  expect_error(.get_formula_vars(abc ~ abc, data = data),
               regexp = 'empty')
  expect_error(.get_formula_vars(abc ~ abc, 'error'), regexp = 'class')
})

##formula
test_that('Testing formula argument', {
  datax <- data.frame(treat = runif(20),
                      pred = rep(c(TRUE, FALSE), 10))
  expect_error(.get_formula_vars(datax), regexp = 'valid')
  expect_error(.get_formula_vars('abc', datax), regexp = 'valid')
  expect_error(.get_formula_vars(abc ~ abc, datax),
               regexp = 'not found')
  expect_error(.get_formula_vars(treat ~ abc * pred, datax),
               regexp = 'columns')
})

#--testing .process_by()--------------------------------------------------------
#test_that('Testing .process_by()', {
#  datax <- data.frame(treat = runif(20),
#                      pred = rep(c(TRUE, FALSE), 10),
#                      pred_fail = c(rep(c(TRUE, FALSE), 9), c(NA, NA)))
#  treat <- treat
#  by <- 'pred'
#  by2 <- NULL
#  by3 <- c()
#
#})

#--testing match_add_args-------------------------------------------------------
test_that('Testing match_add_args()', {
  arglist <- list(formula = formula(treat ~ sex * age),
                  data = data.frame(age = runif(20),
                                    treat= rep(c(TRUE, FALSE), 10),
                                    sex = rep(c('M', 'M', 'F', 'F', 5))),
                  Hess = TRUE, maxit = 120, softmax = TRUE)

  funlist <- list(nnet::multinom, nnet::nnet.default)
  funlist2 <- nnet::multinom
  funlist3 <- list(nnet::multinom, nnet::nnet.formula)

  expect_no_error(match_add_args(arglist, funlist))
  expect_type(match_add_args(arglist, funlist), type = 'list')
  expect_length(match_add_args(arglist, funlist), 18)
  expect_no_error(match_add_args(arglist2, funlist2))
  expect_length(match_add_args(arglist, funlist2), 7)
  expect_no_error(match_add_args(arglist2, funlist3))
})

#--testing scale_0_to_1---------------------------------------------------------
test_that('Testing scale_0_to_1()', {
  data <- data.frame(binary = rep(c(0, 1), 20),
                     binary_lab = rep(c('M', 'F'), 20),
                     numerics = rnorm(20, 156, 12),
                     numerics_01 = runif(20))

expect_no_error(scale_0_to_1(data$binary_lab))
expect_true(is.factor(scale_0_to_1(data$binary_lab)))
expect_length(scale_0_to_1(data$binary_lab), 40)

expect_no_error(lapply(data, scale_0_to_1))
expect_type(lapply(data, scale_0_to_1), 'list')
expect_length(lapply(data, scale_0_to_1), 4)
expect_no_error(as.data.frame(lapply(data, scale_0_to_1)))
})
