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

##extracting terms
