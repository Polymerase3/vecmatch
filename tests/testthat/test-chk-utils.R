#--testing small utils----------------------------------------------------------
test_that("testing chk-utils", {
  # testing .check_gps_matrix
  gps_matrix1 <- data.frame(a = 1)
  expect_error(.check_gps_matrix(gps_matrix1), regexp = "treatment")

  gps_matrix2 <- data.frame(treatment = 1)
  expect_error(.check_gps_matrix(gps_matrix2), regexp = "itself")

  gps_matrix3 <- data.frame(
    treatment = c("a", "b"),
    d = c(1, 2),
    e = c(1, 2)
  )
  expect_error(.check_gps_matrix(gps_matrix3), regexp = "match the unique")

  gps_matrix4 <- data.frame(
    treatment = c("a", "b"),
    a = c(1, 2),
    b = c(1, NA)
  )
  expect_error(.check_gps_matrix(gps_matrix4), regexp = "NA")

  gps_matrix5 <- data.frame(
    treatment = c("a", "b"),
    a = c(1, 2),
    b = c(1, 2)
  )
  expect_error(.check_gps_matrix(gps_matrix5), regexp = "row-wise")

  # testing .check_integer
  expect_error(.check_integer(c("a", 1), x_name = "vector"),
    regexp = "integer"
  )

  # testing match_discrete_args
  expect_error(
    .match_discrete_args("abc",
      choices = c(1, 2, 3),
      x_name = "arg"
    ),
    regexp = "following values"
  )
})
