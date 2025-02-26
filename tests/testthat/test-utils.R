#--testing get_formula_vars()---------------------------------------------------
## data
test_that("Testing data argument", {
  data <- data.frame()
  datax <- data.frame(
    treat = runif(20),
    pred = rep(c(TRUE, FALSE), 10)
  )

  expect_no_error(.get_formula_vars(datax$treat ~ datax$pred))
  expect_no_error(.get_formula_vars(treat ~ pred, data = datax))
  expect_error(.get_formula_vars(abc ~ abc, data = data),
    regexp = "empty"
  )
  expect_error(.get_formula_vars(abc ~ abc, "error"), regexp = "class")
})

## formula
test_that("Testing formula argument", {
  datax <- data.frame(
    treat = runif(20),
    pred = rep(c(TRUE, FALSE), 10)
  )

  expect_error(.get_formula_vars(datax), regexp = "valid")
  expect_error(.get_formula_vars("abc", datax), regexp = "valid")
  expect_error(.get_formula_vars(abc ~ abc, datax),
    regexp = "not found"
  )
  expect_error(.get_formula_vars(treat ~ abc * pred, datax),
    regexp = "columns"
  )

  ## rhs is data.frame
  expect_no_error(.get_formula_vars(datax$treat ~ datax, datax))
})

#--testing .process_by()--------------------------------------------------------


#--testing match_add_args-------------------------------------------------------
test_that("Testing match_add_args()", {
  arglist <- list(
    formula = formula(treat ~ sex * age),
    data = data.frame(
      age = runif(20),
      treat = rep(c(TRUE, FALSE), 10),
      sex = rep(c("M", "M", "F", "F", 5))
    ),
    Hess = TRUE, maxit = 120, softmax = TRUE
  )

  funlist <- list(nnet::multinom, nnet::nnet.default)
  funlist2 <- nnet::multinom
  funlist3 <- list(nnet::multinom, nnet::nnet.formula)

  expect_no_error(match_add_args(arglist, funlist))
  expect_type(match_add_args(arglist, funlist), type = "list")
  expect_length(match_add_args(arglist, funlist), 19)
  expect_no_error(match_add_args(arglist, funlist2))
  expect_length(match_add_args(arglist, funlist2), 7)
  expect_no_error(match_add_args(arglist, funlist3))
})

#--testing scale_0_to_1---------------------------------------------------------
test_that("Testing scale_0_to_1()", {
  data <- data.frame(
    binary = rep(c(0, 1), 20),
    binary_lab = rep(c("M", "F"), 20),
    numerics = rnorm(20, 156, 12),
    numerics_01 = runif(20)
  )

  expect_no_error(scale_0_to_1(data$binary_lab))
  expect_true(is.factor(scale_0_to_1(data$binary_lab)))
  expect_length(scale_0_to_1(data$binary_lab), 40)

  expect_no_error(lapply(data, scale_0_to_1))
  expect_type(lapply(data, scale_0_to_1), "list")
  expect_length(lapply(data, scale_0_to_1), 4)
  expect_no_error(as.data.frame(lapply(data, scale_0_to_1)))
})

#--testing other small utils----------------------------------------------------
test_that("Testing small utils", {
  # make_list
  expect_error(.make_list(NULL), regexp = "integer")

  # nunique
  expect_equal(nunique(NULL), 0)
  expect_equal(nunique(c("a", NA)), 1)

  # na.rem
  expect_no_error(na_rem(c(1:5, NA)))

  # word_list
  expect_equal(attr(word_list(), "plural"), FALSE)

  compare1 <- "a is"
  attr(compare1, "plural") <- FALSE
  expect_equal(word_list("a", is_are = TRUE), compare1)

  compare2 <- "a, b, c are"
  attr(compare2, "plural") <- TRUE
  expect_equal(
    word_list(c("a", "b", "c"), is_are = TRUE, and_or = NULL),
    compare2
  )

  # add_quotes
  expect_equal(add_quotes("a", quotes = TRUE), "\"a\"")
  expect_equal(add_quotes("a", quotes = "b"), "bab")
  expect_error(add_quotes("a", quotes = c("a", "b")), regexp = "boolean")
  expect_equal(add_quotes("a", quotes = 0), "a")
  expect_equal(add_quotes("a", quotes = 1), "'a'")

  # .process_by
  data <- data.frame(
    treat = rep(c("A", "B", "C"), 7),
    y = runif(21),
    sex = c(rep(c("M", "F", "F", "M"), 5), "M"),
    sex_fail = c(rep(c("M", "F", "F", "M"), 5), NA),
    sex_fail2 = c(rep(c("M", "M", "M", "M"), 5), "F")
  )

  expect_error(.process_by(data = data, treat = "treat"),
    regexp = "single string"
  )
  expect_error(.process_by(by = "sex_fail", data = data, treat = "treat"),
    regexp = "contain any"
  )
  expect_error(.process_by(by = "sex_fail2", data, "treat"),
    regexp = "all treatment levels"
  )

  # all_the_same
  expect_equal(all_the_same(c("x", NA), na_rm = FALSE), FALSE)
  expect_equal(all_the_same(c("x", NA), na_rm = TRUE), TRUE)

  # scale_0_to_1
  x <- c(TRUE, FALSE)
  y <- factor(c("a", "a"))
  expect_equal(scale_0_to_1(x), x)
  expect_equal(scale_0_to_1(y), y)
  expect_error(
    scale_0_to_1(c(
      as.Date("01.01.2025", "%d.%m.%Y"),
      as.Date("02.01.2025", "%d.%m.%Y"),
      as.Date("03.01.2025", "%d.%m.%Y")
    )),
    regexp = "Cannot convert"
  )

  # %nin%
  expect_equal(1 %nin% 0, TRUE)
})
