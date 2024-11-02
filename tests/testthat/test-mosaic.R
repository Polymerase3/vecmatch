# test
test_that("multiplication works", {
  # Generate a synthetic dataset
  data <- data.frame(
    age = sample(c("18-25", "26-35", "36-45"), 100, replace = TRUE),
    sex = sample(c(0, 1), 100, replace = TRUE),
    product = sample(c("Electronics", "Clothing", "Food"), 100, replace = TRUE)
  )
})
