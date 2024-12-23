# --fast checking formals-------------------------------------------------------
test_that("mosaic plot: checking formals", {
  withr::with_seed("164678", {
    # Generate a synthetic dataset
    data <- data.frame(
      age = sample(c("18-25", "26-35", "36-45"), 1000,
        replace = TRUE
      ),
      sex = sample(c(0, 1), 1000, replace = TRUE),
      product = sample(c("Electronics", "Clothing", "Food"),
        1000,
        replace = TRUE
      )
    )
  })
  data_fail <- data.frame(y = c("a", "b", "c"))

  # check data and mapping
  expect_error(mosaic(data_fail), regexp = "missing")
  expect_error(mosaic(data_fail, y, fail), regexp = "not in")

  # check single aes run
  expect_no_error(mosaic(data, age))
  expect_error(mosaic(data, age, significance = TRUE), regexp = "group")
  expect_error(
    mosaic(data[data$age == "18-25", ], sex, age,
      significance = TRUE
    ),
    regexp = "impossible"
  )

  # test plot_name arg
  withr::with_tempdir({
    mosaic(data, age, sex, plot_name = "myplot.png")

    # expect that the file exists
    expect_true(file.exists("myplot.png"))

    # expect error when not allowed to overwrite
    expect_error(mosaic(data, age, sex, plot_name = "myplot.png"),
      regexp = "exists"
    )

    # expect no error when allowed to overwrite
    expect_no_error(mosaic(data, age, sex,
      plot_name = "myplot.png",
      overwrite = TRUE
    ))
  })
})

# --checking statistical options------------------------------------------------
test_that("mosai plots: statistics", {
  withr::with_seed("164678", {
    # Generate a synthetic dataset
    data <- data.frame(
      age = sample(c("18-25", "26-35", "36-45"),
        1000,
        replace = TRUE
      ),
      sex = sample(c(0, 1), 1000, replace = TRUE),
      product = sample(c("Electronics", "Clothing", "Food"),
        1000,
        replace = TRUE
      )
    )
  })

  ## check group counts
  expect_no_error(mosaic(data, age, sex, product,
    significance = TRUE,
    group_counts = TRUE
  ))
  expect_no_error(mosaic(data, age, sex,
    significance = TRUE,
    group_counts = TRUE
  ))
  expect_no_error(mosaic(data, age, facet = sex))

  ## to not count as empty test
  expect_error(mosaic())
})
