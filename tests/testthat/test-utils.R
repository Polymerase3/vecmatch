#--testing get_formula_vars()---------------------------------------------------
##data
test_that('Testing data argument', {
  data <- data.frame()
  datax <- data.frame(random = double())
  expect_no_error(.get_formula_vars(abc ~ abc))
  expect_no_error(.get_formula_vars(abc ~ abc, data = datax))
  expect_error(.get_formula_vars(abc ~ abc, data = data),
               regexp = 'empty')
  expect_error(.get_formula_vars(abc ~ abc, 'error'), regexp = 'class')
})

##formula
test_that('Testing formula argument', {
  datax <- data.frame(random = double())
  expect_error(.get_formula_vars(datax), regexp = 'valid')
  expect_error(.get_formula_vars('abc', datax), regexp = 'valid')
  expect_no_error(.get_formula_vars(abc ~ abc, datax))
  expect_no_error(.get_formula_vars(abc ~ abc))

})

##extracting terms
