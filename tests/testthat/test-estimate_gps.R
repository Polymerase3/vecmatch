## --testing formals: formula---------------------------------------------------
test_that('Formals checking: formula', {
  expect_error(estimate_gps(), regexp = 'missing')
  expect_error(estimate_gps(NULL), regexp = 'treatment')
})














