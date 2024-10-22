#-- checking estimate_gps_multinom----------------------------------------------
#-- 'multinom' already entirely checked-----------------------------------------

#-- checking vglm (VGAM package)------------------------------------------------
test_that("Checking `vglm` method", {
  data <- data.frame(treat = rep(c(1, 2, 3, 4, 5), 20),
                     pred = runif(100),
                     pred2 = runif(100))

  # link
  expect_error(estimate_gps(treat ~ pred, data, method = 'vglm', link = 'multinom'),
               regexp = 'only accepts')
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'vglm', link = 'multinomial_logit'))
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'vglm', link = 'reduced_rank_ml'))

  # control processing
  expect_error(estimate_gps(treat ~ pred, data, method = 'vglm',
                            control = 'failure', maker = 123), regexp = 'valid')
  expect_error(estimate_gps(treat ~ pred, data, method = 'vglm',
                            control = 'failure', link = 'reduced_rank_ml',
                            maker = 123), regexp = 'valid')
  expect_no_error(estimate_gps(treat ~ pred * pred2, data, method = 'vglm',
                            control = VGAM::vglm.control(), maker = 'pass'))
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'vglm',
                               control = VGAM::vglm.control(step = 2, maxit = 50),
                               maker = 'pass'))

  # family processing
  expect_error(estimate_gps(treat ~ pred, data, method = 'vglm', family = 'fail'),
               regexp = 'VGAM family')
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'vglm', family = VGAM::multinomial()))

})

##-- checking brglm2 package---------------------------------------------------
test_that('Checking `brglm2` method', {
  data <- data.frame(treat = rep(c(1, 2, 3, 4, 5), 20),
                     pred = runif(100),
                     pred2 = runif(100))

  # basic run
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'brglm2'))

  # testing control
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'brglm2',
                               control = brglm2::brglmControl()))

  # testing other args in control
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'brglm2',
                               control = brglm2::brglmControl(slowit = 1, check_aliasing = FALSE)))
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'brglm2',
                               control = brglm2::brglmControl(slowit = 1, check_aliasing = FALSE)))
  expect_error(estimate_gps(treat ~ pred, data, method = 'brglm2', model = TRUE,
                               control = 'fail'), regexp = 'valid function call')


})

##-- checking brglm2 package---------------------------------------------------
test_that('Checking `mclogit` method', {
  data <- data.frame(treat = rep(c(1, 2, 3, 4, 5), 20),
                     pred = runif(100),
                     pred2 = runif(100),
                     sex = rep(c(1, 2), 50))

  # basic run
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'mblogit'))

  # link fun
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'mblogit',
                               link = 'baseline_category_logit'))
  expect_error(estimate_gps(treat ~ pred, data, method = 'mblogit', link = 'fail'),
               regexp = 'only accepts')

  # control testing
  expect_no_error(estimate_gps(treat ~ pred, data, method = 'mblogit',
                               control = mclogit::mclogit.control(maxit = 30)))
  expect_error(estimate_gps(treat ~ pred, data, method = 'mblogit',
                            control = 'fail'))
})
