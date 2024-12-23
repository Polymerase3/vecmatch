test_that("testing plot utils", {
  # .generate_colors
  output_list <- list(
    c("black"),
    c("black", "red"),
    c("black", "red", "white")
  )
  expect_equal(
    .generate_colors(c("black", "red", "white")),
    output_list
  )

  # scales
  expect_contains(
    class(scale_fill_vecmatch(1, type = "continuous")),
    "ScaleContinuous"
  )
  expect_contains(
    class(scale_color_vecmatch(1, type = "continuous")),
    "ScaleContinuous"
  )
})
