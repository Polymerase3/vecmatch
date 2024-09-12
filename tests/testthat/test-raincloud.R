## --testing formals: data-------------------------------------------------
test_that('Formals checking: data', {
  expect_error(raincloud(), regexp = 'data.frame')

  data <- character()
  expect_error(raincloud(data), regexp = 'class')

  data <- data.frame()
  expect_error(raincloud(data), regexp = 'empty')

  data <- data.frame(x = character())
  expect_error(raincloud(data), regexp = 'numeric', label = 'not numeric column')
})

## --testing formals: y, group, facet-------------------------------------------
test_that('Formals checking: data', {
  data <- data.frame(x = double())
  expect_error(raincloud(data, c(1, 2)), regexp = 'column names')
  expect_no_error(raincloud(data, x))
  expect_no_error(raincloud(data, 'x'))
})

## --testing if names in the names(data)----------------------------------------
## --testitng .check_name
test_that('Formals checking: names provided in the colnames', {
  data <- data.frame(name1 = numeric(),
                     name2 = numeric(),
                     name3 = factor(),
                     string = character())
  fail <- list('name1', 'name2', 'fail', 'fail')
  success <- list('name1', 'name2', 'name3', 'string')

  expect_equal(.check_name(data, fail), c('fail', 'fail'))
  expect_invisible(.check_name(data, success))

  raincloud(data, y = )


})
