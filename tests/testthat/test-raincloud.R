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
  data <- data.frame(random = double())
  expect_error(raincloud(data, c(1, 2)), regexp = 'valid')
  expect_error(raincloud(data, facet = random), regexp = 'default')
  expect_no_error(raincloud(data, y = random))
  expect_no_error(raincloud(data, y = 'random'))
})

## --testing if names in the names(data)----------------------------------------
## --testitng .check_name
test_that('Formals checking: names provided in the colnames', {
  data <- data.frame(name1 = numeric(),
                     name2 = numeric(),
                     name3 = factor(),
                     string = character())
  expect_error(raincloud(data), regexp = 'default')

  fail <- list('name1', 'name2', 'fail', 'fail')
  success <- list('name1', 'name2', 'name3', 'string')

  expect_equal(.check_name(data, fail), c('fail', 'fail'))
  expect_no_error(.check_name(data, success))

  expect_error(raincloud(data, y = random), regexp = 'random')
  expect_error(raincloud(data, y = 'random'), regexp = 'random')
  expect_error(raincloud(data, y = 'name1', facet = 'string', group = random),
               regexp = 'random')
  expect_no_error(raincloud(data, y = 'name1', facet = name2, group = string))
})

## --testing significance-------------------------------------------------------
test_that('Formals checking: significance', {
  data <- data.frame(random = double())
  expect_error(raincloud(data, random, significance = 'asd'),
                regexp = 'logical')
})

## --testing limits-------------------------------------------------------------
test_that('Formals checking: limits', {
  data <- data.frame(random = double())
  expect_error(raincloud(data, random, limits = 'asd'),
               regexp = 'limits')
  expect_error(raincloud(data, random, limits = c(1, 2, 3)),
               regexp = 'limits')
  expect_error(raincloud(data, random, limits = c(1, 'asd')),
               regexp = 'limits')
  expect_no_error(raincloud(data, random, limits = c(1, 2)))
})

## --testing jitter-------------------------------------------------------------
test_that('Formals checking: jitter', {
  data <- data.frame(random = double())
  expect_error(raincloud(data, random, jitter = 2), regexp = 'between')
  expect_no_error(raincloud(data, random, jitter = 0.7))
})

## -testing alpha---------------------------------------------------------------
test_that('Formals checking: alpha', {
  data <- data.frame(random = double())
  expect_error(raincloud(data, random, alpha = 2), regexp = 'between')
  expect_no_error(raincloud(data, random, alpha = 0.7))
})

## --testing significance-------------------------------------------------------
test_that('Formals checking: save', {
  data <- data.frame(random = double())
  expect_error(raincloud(data, random, save = 'asd'),
               regexp = 'logical')
})

## --testing non-numeric y------------------------------------------------------
test_that('Data converting: numeric', {
  data <- data.frame(pass = c(1, 2, 3),
                     fail = c('a', 'b', 'c'),
                     fail2 = c('1', 'TRUE', 'c'),
                     pass2 = c('1', '2', '3'),
                     pass3 = c(TRUE, FALSE, FALSE))
  expect_error(raincloud(data, y = fail), regexp = 'numeric')
  expect_error(raincloud(data, y = fail2), regexp = 'numeric')
  expect_no_error(raincloud(data, y = 'pass'))
  expect_no_error(raincloud(data, y = pass2))
  expect_no_error(raincloud(data, y = 'pass3'))
})

## --testing factor conversion--------------------------------------------------
test_that('Data converting: factors', {
  data <- data.frame(y = runif(20),
                     charr = rep(c('a', 'b', 'c', 'd'), 5),
                     logg = rep(c(TRUE, FALSE), 10),
                     int = rep(c(1, 2, 3, 4, 5), 4),
                     fact = factor(rep(c(1, 2, 3, 4), 5)),
                     fail1 = c(rep(as.Date('01-01-2021'), 20)),
                     warning = 1:20)

  expect_no_error(raincloud(data, y, facet = charr))
  expect_no_error(raincloud(data, y = y, group = 'charr'))
  expect_no_error(raincloud(data, y = y, group = charr, facet = 'logg'))
  expect_no_error(raincloud(data, y, facet = logg))
  expect_no_error(raincloud(data, y = y, group = logg))
  expect_no_error(raincloud(data, y = 'y', group = logg, facet = int))
  expect_no_error(raincloud(data, y, facet = 'int'))
  expect_no_error(raincloud(data, y = y, group = int))
  expect_no_error(raincloud(data, y = y, group = 'int', facet = 'fact'))
  expect_error(raincloud(data, y, fail1), regexp = 'converted')
  expect_warning(raincloud(data, y, facet = 'warning'), regexp = '10')
})
