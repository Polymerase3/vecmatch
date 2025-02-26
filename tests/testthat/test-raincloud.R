## --testing formals: datax-------------------------------------------------
test_that("Formals checking: data", {
  expect_error(raincloud(), regexp = "data.frame")

  datax <- character()
  expect_error(raincloud(datax), regexp = "class")

  datax <- data.frame()
  expect_error(raincloud(datax), regexp = "empty")

  datax <- data.frame(x = character())
  expect_error(raincloud(datax), regexp = "numeric")

  datax <- data.frame(y = runif(100))
  expect_no_error(raincloud(datax, y))
})

## --testing formals: y, group, facet-------------------------------------------
test_that("Formals checking: data", {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, c(1, 2)), regexp = "valid")
  expect_error(raincloud(datax, facet = random), regexp = "default")
  expect_no_error(raincloud(datax, y = random))
  expect_no_error(raincloud(datax, y = "random"))
})

## --testing if names in the names(datax)---------------------------------------
## --testitng .check_name
test_that("Formals checking: names provided in the colnames", {
  datax <- data.frame(
    name1 = numeric(),
    name2 = numeric(),
    name3 = factor(),
    string = character()
  )

  expect_error(raincloud(datax), regexp = "default")

  fail <- list("name1", "name2", "fail", "fail")
  success <- list("name1", "name2", "name3", "string")

  expect_equal(.check_name(datax, fail), c("fail", "fail"))
  expect_no_error(.check_name(datax, success))

  expect_error(raincloud(datax, y = "name1", facet = "string", group = random),
    regexp = "random"
  )
  expect_error(raincloud(datax, y = "name1", facet = name2, group = string),
    regexp = "observations"
  )
})

## --testing significance-------------------------------------------------------
test_that("Formals checking: significance", {
  datax <- data.frame(random = double())

  datax2 <- data.frame(
    y = runif(10),
    charr = rep(c("a", "b", "c", "d"), 25),
    logg = sample(c(0, 1), size = 100, replace = TRUE)
  )

  expect_error(raincloud(datax, random, significance = "asd"),
    regexp = "one group"
  )
  raincloud(datax2, y, logg, significance = "t_test")
  expect_no_error(raincloud(datax2, y, logg, significance = "t_test"))
  expect_no_error(raincloud(datax2, y, logg, significance = "tukeyHSD_test"))
  expect_error(raincloud(datax2, y, logg, significance = "test"),
    regexp = "methods"
  )
  expect_no_error(raincloud(datax2, y, charr, logg, significance = "t_test"))
  expect_no_error(raincloud(datax2, y, charr, logg,
    significance = "t_test",
    limits = c(-1, 2)
  ))
})

## --testing limits-------------------------------------------------------------
test_that("Formals checking: limits", {
  datax <- data.frame(random = double())
  datax2 <- data.frame(y = runif(20))

  expect_error(raincloud(datax, random, limits = "asd"),
    regexp = "limits"
  )
  expect_error(raincloud(datax, random, limits = c(1, 2, 3)),
    regexp = "limits"
  )
  expect_error(raincloud(datax, random, limits = c(1, "asd")),
    regexp = "limits"
  )
  expect_no_error(raincloud(datax2, y, limits = c(1, 2)))
})

## --testing jitter-------------------------------------------------------------
test_that("Formals checking: jitter", {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, random, jitter = 2), regexp = "between")
  expect_no_error(raincloud(datax, random, jitter = 0.7))
})

## -testing alpha---------------------------------------------------------------
test_that("Formals checking: alpha", {
  datax <- data.frame(random = double())
  expect_error(raincloud(datax, random, alpha = 2), regexp = "between")
  expect_no_error(raincloud(datax, random, alpha = 0.7))
})

## --testing plot_name----------------------------------------------------------
test_that("Formals checking: plot_name", {
  datax <- data.frame(
    y = runif(20),
    group = rep(c(TRUE, FALSE), 10)
  )

  expect_error(raincloud(datax, y, plot_name = c(1, 2)),
    regexp = "character"
  )
  expect_error(raincloud(datax, y, plot_name = "invalid"),
    regexp = ".pdf"
  )
  expect_error(raincloud(datax, y, plot_name = "invalid.sav"),
    regexp = ".pdf"
  )

  withr::with_tempdir({
    expect_no_error(raincloud(datax, y,
      plot_name = "valid.png",
      overwrite = TRUE
    ))
  })
})

## --testing overwrite and sig_label_color--------------------------------------
test_that("Formals checking: overwrite and sig_label_color", {
  datax <- data.frame(
    y = runif(20),
    group = rep(c(TRUE, FALSE), 10)
  )

  expect_error(raincloud(datax, y, overwrite = "asd"),
    regexp = "flag"
  )

  expect_no_error(raincloud(datax, y, group,
    significance = "t_test",
    sig_label_color = TRUE
  ))
})

## --testing smd_type-----------------------------------------------------------
test_that("Formals checking: smd_type", {
  datax <- data.frame(
    y = runif(20),
    group = rep(c(TRUE, FALSE), 10)
  )

  expect_error(raincloud(datax, y, smd_type = NULL), regexp = "character")
  expect_error(raincloud(datax, y, smd_type = c("a", "b")), regexp = "length")
  expect_error(raincloud(datax, y, smd_type = "sd"), regexp = "mean")

  expect_no_error(raincloud(datax, y, group, smd_type = "median"))
})

## --testing non-numeric y------------------------------------------------------
test_that("datax converting: numeric", {
  datax <- data.frame(
    pass = c(1, 2, 3),
    fail = c("a", "b", "c"),
    pass2 = c("1", "2", "3"),
    pass3 = c(TRUE, FALSE, FALSE)
  )

  expect_error(raincloud(datax, y = fail), regexp = "numeric")
  expect_no_error(raincloud(datax, y = "pass"))
  expect_no_error(raincloud(datax, y = pass2))
  expect_no_error(raincloud(datax, y = "pass3"))
})

## --testing factor conversion--------------------------------------------------
test_that("datax converting: factors", {
  datax <- data.frame(
    y = runif(20),
    charr = rep(c("a", "b", "c", "d"), 5),
    logg = rep(c(TRUE, FALSE), 10),
    int = rep(c(1, 2, 3, 4, 5), 4),
    fact = factor(rep(c(1, 2, 3, 4), 5)),
    fail1 = c(rep(as.Date("01-01-2021"), 20)),
    warning = 1:20
  )

  expect_no_error(raincloud(datax, y, facet = charr))
  expect_no_error(raincloud(datax, y = y, group = "charr"))
  expect_no_error(raincloud(datax, y = y, group = charr, facet = "logg"))
  expect_no_error(raincloud(datax, y, facet = logg))
  expect_no_error(raincloud(datax, y = y, group = logg))
  expect_no_error(raincloud(datax, y = "y", group = logg, facet = int))
  expect_no_error(raincloud(datax, y, facet = "int"))
  expect_no_error(raincloud(datax, y = y, group = int))
  expect_no_error(raincloud(datax, y = y, group = "int", facet = "fact"))
  expect_error(raincloud(datax, y, fail1), regexp = "converted")
  expect_warning(raincloud(datax, y, facet = "warning"), regexp = "10")
})

## --testing significance methods-----------------------------------------------
test_that("significance: testing methods", {
  datax <- data.frame(
    y = runif(100),
    charr = rep(c("a", "b"), 50)
  )

  expect_no_error(raincloud(datax, y, charr, significance = "t_test"))
  expect_no_error(raincloud(datax, y, charr, significance = "wilcoxon_test"))
  expect_no_error(raincloud(datax, y, charr, significance = "dunn_test"))
  expect_no_error(raincloud(datax, y, charr, significance = "tukeyHSD_test"))
  expect_no_error(raincloud(datax, y, charr,
    significance = "games_howell_test"
  ))
  expect_no_error(raincloud(datax, y, charr, significance = "sign_test"))
  expect_error(raincloud(datax, y, charr, significance = "some_method"),
    regexp = "methods"
  )
})

## --testing sig_label_size-----------------------------------------------------
test_that("sig_label_size: argument check", {
  datax <- data.frame(y = runif(20))

  expect_no_error(raincloud(datax, y, sig_label_size = 8))
  expect_no_error(raincloud(datax, y, sig_label_size = 8L))
  expect_error(raincloud(datax, y, sig_label_size = "a"), regexp = "integer")
})

## --testing plot_name and overwrite--------------------------------------------
test_that("plot_name and overwrite: argument and functioning check", {
  datax <- data.frame(
    y = runif(100),
    charr = rep(c("a", "b"), 50)
  )

  withr::with_tempdir({
    raincloud(datax, y, charr, plot_name = "myplot.png")

    # expect that the file exists
    expect_true(file.exists("myplot.png"))

    # expect error when not allowed to overwrite
    expect_error(raincloud(datax, y, charr, plot_name = "myplot.png"),
      regexp = "exists"
    )



    # expect no error when allowed to overwrite
    expect_no_error(raincloud(datax, y, charr,
      plot_name = "myplot.png",
      overwrite = TRUE
    ))
  })

  # control if cleaned
  expect_false(file.exists("myplot.png"))
})
