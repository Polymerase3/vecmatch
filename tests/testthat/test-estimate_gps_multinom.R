#-- checking estimate_gps_multinom----------------------------------------------
#-- 'multinom' already entirely checked-----------------------------------------

#-- checking vglm (VGAM package)------------------------------------------------
test_that("Checking `vglm` method", {
  data <- data.frame(treat = rep(c(1, 2, 3, 4, 5), 20),
                     pred = runif(100),
                     pred2 = runif(100))


})
