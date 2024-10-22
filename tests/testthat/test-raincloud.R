## --testing formals: datax-------------------------------------------------
test_that('Formals checking: data', {
  expect_error(raincloud(), regexp = 'data.frame')

  datax <- character()
  expect_error(raincloud(datax), regexp = 'class')

  datax <- data.frame()
  expect_error(raincloud(datax), regexp = 'empty')

  datax <- data.frame(x = character())
  expect_error(raincloud(datax), regexp = 'numeric', label = 'not numeric column')
})

## --testing formals: y, group, facet-------------------------------------------
test_that('Formals checking: data', {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, c(1, 2)), regexp = 'valid')
  expect_error(raincloud(datax, facet = random), regexp = 'default')
  expect_no_error(raincloud(datax, y = random))
  expect_no_error(raincloud(datax, y = 'random'))
})

## --testing if names in the names(datax)----------------------------------------
## --testitng .check_name
test_that('Formals checking: names provided in the colnames', {
  datax <- data.frame(name1 = numeric(),
                     name2 = numeric(),
                     name3 = factor(),
                     string = character())
  expect_error(raincloud(datax), regexp = 'default')

  fail <- list('name1', 'name2', 'fail', 'fail')
  success <- list('name1', 'name2', 'name3', 'string')

  expect_equal(.check_name(datax, fail), c('fail', 'fail'))
  expect_no_error(.check_name(datax, success))

  expect_error(raincloud(datax, y = random), regexp = 'random')
  expect_error(raincloud(datax, y = 'random'), regexp = 'random')
  expect_error(raincloud(datax, y = 'name1', facet = 'string', group = random),
               regexp = 'random')
  expect_no_error(raincloud(datax, y = 'name1', facet = name2, group = string))
})

## --testing significance-------------------------------------------------------
test_that('Formals checking: significance', {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, random, significance = 'asd'),
                regexp = 'logical')
})

## --testing limits-------------------------------------------------------------
test_that('Formals checking: limits', {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, random, limits = 'asd'),
               regexp = 'limits')
  expect_error(raincloud(datax, random, limits = c(1, 2, 3)),
               regexp = 'limits')
  expect_error(raincloud(datax, random, limits = c(1, 'asd')),
               regexp = 'limits')
  expect_no_error(raincloud(datax, random, limits = c(1, 2)))
})

## --testing jitter-------------------------------------------------------------
test_that('Formals checking: jitter', {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, random, jitter = 2), regexp = 'between')
  expect_no_error(raincloud(datax, random, jitter = 0.7))
})

## -testing alpha---------------------------------------------------------------
test_that('Formals checking: alpha', {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, random, alpha = 2), regexp = 'between')
  expect_no_error(raincloud(datax, random, alpha = 0.7))
})

## --testing save---------------------------------------------------------------
test_that('Formals checking: save', {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, random, save = 'asd'),
               regexp = 'logical')
})

## --testing plot.name----------------------------------------------------------
test_that('Formals checking: plot.name', {
  datax <- data.frame(y = runif(20),
                      group = rep(c(TRUE, FALSE), 10))
  expect_error(raincloud(datax, y, save = TRUE, plot.name = c(1, 2)),
               regexp = 'character')
  expect_error(raincloud(datax, y, save = TRUE),
               regexp = 'specified')
  expect_error(raincloud(datax, y, save = TRUE,
                         plot.name = 'invalid'),
               regexp = '.pdf')
  expect_error(raincloud(datax, y, save = TRUE,
                         plot.name = 'invalid.sav'),
               regexp = '.pdf')
  expect_no_error(raincloud(datax, y, save = TRUE,
                            plot.name = 'valid.png'))
  if(file.exists('valid.png')) file.remove('valid.png')
})
## --testing overwrite----------------------------------------------------------
test_that('Formals checking: overwrite', {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, random, overwrite = 'asd'),
               regexp = 'logical')
})

## --testing non-numeric y------------------------------------------------------
test_that('datax converting: numeric', {
  datax <- data.frame(pass = c(1, 2, 3),
                     fail = c('a', 'b', 'c'),
                     fail2 = c('1', 'TRUE', 'c'),
                     pass2 = c('1', '2', '3'),
                     pass3 = c(TRUE, FALSE, FALSE))
  expect_error(raincloud(datax, y = fail), regexp = 'numeric')
  expect_error(raincloud(datax, y = fail2), regexp = 'numeric')
  expect_no_error(raincloud(datax, y = 'pass'))
  expect_no_error(raincloud(datax, y = pass2))
  expect_no_error(raincloud(datax, y = 'pass3'))
})

## --testing factor conversion--------------------------------------------------
test_that('datax converting: factors', {
  datax <- data.frame(y = runif(20),
                     charr = rep(c('a', 'b', 'c', 'd'), 5),
                     logg = rep(c(TRUE, FALSE), 10),
                     int = rep(c(1, 2, 3, 4, 5), 4),
                     fact = factor(rep(c(1, 2, 3, 4), 5)),
                     fail1 = c(rep(as.Date('01-01-2021'), 20)),
                     warning = 1:20)

  expect_no_error(raincloud(datax, y, facet = charr))
  expect_no_error(raincloud(datax, y = y, group = 'charr'))
  expect_no_error(raincloud(datax, y = y, group = charr, facet = 'logg'))
  expect_no_error(raincloud(datax, y, facet = logg))
  expect_no_error(raincloud(datax, y = y, group = logg))
  expect_no_error(raincloud(datax, y = 'y', group = logg, facet = int))
  expect_no_error(raincloud(datax, y, facet = 'int'))
  expect_no_error(raincloud(datax, y = y, group = int))
  expect_no_error(raincloud(datax, y = y, group = 'int', facet = 'fact'))
  expect_error(raincloud(datax, y, fail1), regexp = 'converted')
  expect_warning(raincloud(datax, y, facet = 'warning'), regexp = '10')
})

## --testing significance methods-----------------------------------------------
