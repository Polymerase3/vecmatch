test_that("Testing make_opt_args S3 class", {
  data <- data.frame(
    treat = rep(c(1, 2, 3, 4), 25),
    y     = rep(c(TRUE, FALSE), 10),
    pred  = rnorm(20)
  )

  formula_test <- formula(treat ~ pred * y)

  # check data
  ## not a data.frame
  expect_error(make_opt_args(list()), "data.frame")

  # check formula
  ## no formula
  expect_error(make_opt_args(), "formula")

  ## rest is standard workflow already checked in estimate_gps()

  # check gps_methods
  expect_error(
    make_opt_args(
      data,
      formula_test,
      gps_method = "abc"
    ),
    "m1"
  )

  expect_error(
    make_opt_args(data, formula_test, gps_method = c("m1", "m1")),
    "Duplicates"
  )

  expect_error(
    make_opt_args(
      data,
      formula_test,
      gps_method = paste0("m", 1:11)
    ),
    "length"
  )

  # check reference
  ## reference not in treatment levels
  expect_error(
    make_opt_args(data, formula_test, "m1", reference = "BC"),
    "unique levels"
  )

  ## reference not a string
  expect_error(
    make_opt_args(data, formula_test, "m1", reference = 3),
    "string"
  )

  # matching method
  expect_error(
    make_opt_args(data, formula_test, matching_method = "abc"),
    "matching_method"
  )

  # caliper
  ## caliper is 0
  expect_error(
    make_opt_args(data, formula_test, caliper = 0),
    "Zeros"
  )

  ## caliper is not numeric
  expect_error(
    make_opt_args(data, formula_test, caliper = "a"),
    "numeric"
  )

  ## caliper duplicated
  expect_error(
    make_opt_args(data, formula_test, caliper = c(1, 1, 1.2)),
    "Duplicates"
  )

  # order
  expect_error(
    make_opt_args(data, formula_test, order = "ab"),
    "order"
  )

  # cluster
  ## cluster exceeds number of unique treatment levels
  expect_error(make_opt_args(data, formula_test, cluster = 5), "cluster")

  # method = "nnm"
  ## validate ratio
  expect_error(
    make_opt_args(
      data,
      formula_test,
      matching_method = "nnm", ratio = 21
    ),
    "ratio"
  )

  ## validate replace
  expect_error(
    make_opt_args(
      data,
      formula_test,
      matching_method = "nnm", replace = "21"
    ),
    "replace"
  )

  ## validate ties
  expect_error(
    make_opt_args(
      data,
      formula_test,
      matching_method = "nnm", ties = "21"
    ),
    "ties"
  )

  # method = "fullopt"
  ## validate min and max controls
  expect_error(
    make_opt_args(
      data,
      formula_test,
      matching_method = "fullopt",
      min_controls    = "21"
    ),
    "min_controls"
  )

  expect_error(
    make_opt_args(
      data,
      formula_test,
      matching_method = "fullopt",
      max_controls    = "21"
    ),
    "max_controls"
  )

  # check manually number of combinations
  ## both methods
  both <- make_opt_args(
    data,
    formula_test,
    reference       = c("1", "3"),
    gps_method      = paste0("m", 1:10),
    matching_method = c("nnm", "fullopt"),
    caliper         = seq(0.37, 2.21, 0.01), # 185
    order           = c("desc", "random", "asc"),
    cluster         = 1:3,
    replace         = c(TRUE, FALSE),
    ties            = FALSE,
    ratio           = 1:3,
    min_controls    = 1:3,
    max_controls    = 3
  )

  expected_nnm <- 2 * 10 * 185 * 3 * 3 * 2 * 3
  expected_fullopt <- 2 * 10 * 185 * 3 * 3 * 3

  expect_equal(
    as.numeric(attr(both, "total_combinations")),
    expected_nnm + expected_fullopt
  )

  ## "nnm"
  nnm <- make_opt_args(
    data,
    formula_test,
    reference       = c("1", "3", "2"),
    gps_method      = c("m1", "m3", "m7", "m10"),
    matching_method = "nnm",
    caliper         = seq(0.1, 2.5, 0.1), # length 25
    order           = c("desc", "random"),
    cluster         = 2:4,
    replace         = TRUE,
    ties            = FALSE,
    ratio           = c(1, 2)
  )

  expected_comb <- 3 * 4 * 25 * 2 * 3 * 2

  expect_equal(as.numeric(attr(nnm, "total_combinations")), expected_comb)

  ## "fullopt"
  fullopt <- make_opt_args(
    data,
    formula_test,
    reference       = c("1", "2", "3", "4"),
    gps_method      = c("m1", "m2", "m3", "m4", "m5", "m6"),
    matching_method = "fullopt",
    caliper         = seq(0.01, 3.45, 0.01), # 345
    order           = "desc",
    cluster         = 1:4,
    min_controls    = 1:2,
    max_controls    = 1:2
  )

  expected_comb <- 4 * 6 * 345 * 4 * 2 * 2

  expect_equal(as.numeric(attr(fullopt, "total_combinations")), expected_comb)

  # print should not error and should show header
  out <- utils::capture.output(print(fullopt))
  expect_true(any(grepl(
    "Optimization Argument Set \\(class: opt_args\\)",
    out
  )))

  ## reproducibility tests -----------------------------------------------------
  # make_opt_args should not change the global rng state
  set.seed(710710)
  old_seed <- .Random.seed
  opt_tmp <- make_opt_args(
    data            = data,
    formula         = formula_test,
    gps_method      = "m1",
    matching_method = "nnm",
    caliper         = c(0.1, 0.2),
    order           = "desc",
    cluster         = 2,
    ratio           = 1,
    replace         = TRUE,
    ties            = TRUE
  )
  expect_identical(.Random.seed, old_seed)

  # make_opt_args run twice with the same seed gives identical object
  set.seed(720720)
  opt1 <- make_opt_args(
    data            = data,
    formula         = formula_test,
    gps_method      = "m1",
    matching_method = "nnm",
    caliper         = c(0.1, 0.2),
    order           = "desc",
    cluster         = 2,
    ratio           = 1,
    replace         = TRUE,
    ties            = TRUE
  )
  set.seed(720720)
  opt2 <- make_opt_args(
    data            = data,
    formula         = formula_test,
    gps_method      = "m1",
    matching_method = "nnm",
    caliper         = c(0.1, 0.2),
    order           = "desc",
    cluster         = 2,
    ratio           = 1,
    replace         = TRUE,
    ties            = TRUE
  )
  expect_identical(opt1, opt2)

  # print.opt_args should not touch rng and be reproducible
  set.seed(730730)
  old_seed <- .Random.seed
  out1 <- utils::capture.output(print(opt1))
  expect_identical(.Random.seed, old_seed)

  set.seed(740740)
  out_a <- utils::capture.output(print(opt1))
  set.seed(740740)
  out_b <- utils::capture.output(print(opt1))
  expect_identical(out_a, out_b)
})

test_that("best_opt_result / select_result methods and helpers work
          end-to-end", {
  # data and formula
  data <- cancer
  formula <- status ~ age * sex

  ## define a relatively small search space to keep runtime reasonable
  opt_args <- make_opt_args(
    data            = data,
    formula         = formula,
    gps_method      = "m1",
    matching_method = "nnm",
    caliper         = c(0.5, 1.0),
    order           = "desc",
    cluster         = 1,
    ratio           = 1,
    replace         = TRUE,
    ties            = TRUE,
    reference       = "control"
  )

  ## run optimization with moderate n_iter
  withr::with_seed(8252, {
    best_res <- optimize_gps(
      data     = data,
      formula  = formula,
      opt_args = opt_args,
      n_iter   = 50
    )
  })

  ## ----- best_opt_result: basic checks --------------------------------------
  expect_s3_class(best_res, "best_opt_result")
  expect_s3_class(best_res, "data.frame")
  expect_true(all(c("smd", "perc_matched") %in% names(best_res)))

  ## summary.best_opt_result ---------------------------------------------------
  s_best <- summary(best_res)
  expect_s3_class(s_best, "summary.best_opt_result")
  expect_s3_class(s_best, "data.frame")
  expect_true(
    all(c("smd_group", "unique_configs", "smd", "perc_matched") %in%
      names(s_best))
  )

  # print.summary.best_opt_result: just check header appears
  out_sum <- utils::capture.output(print(s_best))
  expect_true(
    any(grepl("Best Optimization Results by SMD Group", out_sum, fixed = TRUE))
  )

  ## print.best_opt_result + .print_best_opt_result_core ----------------------
  out_best <- utils::capture.output(
    {
      res_print <- print(best_res)
    },
    type = "message"
  )
  out_best <- cli::ansi_strip(out_best)

  expect_true(
    any(grepl(
      "best_opt_result object \\(GPS matching optimization summary\\)",
      out_best
    ))
  )
  expect_identical(res_print, best_res)

  ## str.best_opt_result -------------------------------------------------------
  out_str <- utils::capture.output(str(best_res))
  expect_true(
    any(grepl(
      "best_opt_result object: GPS matching optimization summary",
      out_str
    ))
  )

  ## plot.best_opt_result ------------------------------------------------------
  expect_silent({
    res_plot <- plot(best_res)
  })
  expect_identical(res_plot, best_res)

  ## ----- select_opt / select_result -----------------------------------------
  # simple selection: use all groups, all covariates, max smd
  sel_res <- select_opt(
    x             = best_res,
    smd_type      = "max",
    smd_groups    = NULL,
    smd_variables = NULL,
    perc_matched  = NULL
  )

  expect_s3_class(sel_res, "select_result")
  expect_s3_class(sel_res, "data.frame")
  expect_true(nrow(sel_res) > 0)
  expect_true(
    all(c("smd_group", "overall_stat", "perc_matched") %in% names(sel_res))
  )

  ## summary.select_result -----------------------------------------------------
  s_sel <- summary(sel_res)
  expect_s3_class(s_sel, "summary.select_result")
  expect_s3_class(s_sel, "data.frame")
  expect_true(
    all(c("smd_group", "unique_configs", "smd", "perc_matched") %in%
      names(s_sel))
  )

  out_s_sel <- utils::capture.output(print(s_sel))
  expect_true(
    any(grepl("Optimization Selection Summary", out_s_sel, fixed = TRUE))
  )

  ## print.select_result -------------------------------------------------------
  out_sel <- utils::capture.output(
    {
      res_sel_print <- print(sel_res)
    },
    type = "message"
  )
  out_sel <- cli::ansi_strip(out_sel)

  expect_true(
    any(grepl(
      "select_result object \\(selected matching configurations\\)",
      out_sel
    ))
  )
  expect_identical(res_sel_print, sel_res)

  ## str.select_result ---------------------------------------------------------
  out_str_sel <- utils::capture.output(str(sel_res))
  expect_true(
    any(grepl(
      "select_result object: selected GPS matching configurations",
      out_str_sel
    ))
  )

  ## get_select_params() -------------------------------------------------------
  params_all <- get_select_params(sel_res)
  expect_s3_class(params_all, "data.frame")
  expect_true(nrow(params_all) > 0)
  expect_true("iter_ID" %in% names(params_all))

  # take one smd_group from selection and filter
  some_group <- as.character(stats::na.omit(sel_res$smd_group)[1L])

  params_group <- get_select_params(sel_res, smd_group = some_group)
  expect_true(nrow(params_group) > 0)
  expect_true(all(params_group$smd_group == some_group))

  ## run_selected_matching() ---------------------------------------------------
  # use the same data/formula as in the optimization, filter by chosen smd_group
  withr::with_seed(123, {
    matched_obj <- run_selected_matching(
      x         = sel_res,
      data      = data,
      formula   = formula,
      smd_group = some_group,
      row       = 1L
    )
  })

  expect_s3_class(matched_obj, "matched")
  expect_s3_class(matched_obj, "data.frame")
  expect_true(!is.null(attr(matched_obj, "original_data")))
  expect_true(!is.null(attr(matched_obj, "matching_filter")))
  expect_true(!is.null(attr(matched_obj, "select_row")))
  expect_true(!is.null(attr(matched_obj, "select_smd_group")))
  expect_true(!is.null(attr(matched_obj, "select_call")))

  ## reproducibility tests -----------------------------------------------------
  # optimize_gps should not change the global rng state and be reproducible
  set.seed(810810)
  old_seed <- .Random.seed
  tmp_opt <- optimize_gps(
    data     = data,
    formula  = formula,
    opt_args = opt_args,
    n_iter   = 5
  )
  expect_identical(.Random.seed, old_seed)

  set.seed(820820)
  opt_a <- optimize_gps(
    data     = data,
    formula  = formula,
    opt_args = opt_args,
    n_iter   = 5
  )
  set.seed(820820)
  opt_b <- optimize_gps(
    data     = data,
    formula  = formula,
    opt_args = opt_args,
    n_iter   = 5
  )

  strip_time <- function(x) {
    attr(x, "optimization_time") <- NULL
    x
  }

  expect_identical(strip_time(opt_a), strip_time(opt_b))

  # summary.best_opt_result should not touch rng and be reproducible
  set.seed(830830)
  old_seed <- .Random.seed
  s_best2 <- summary(best_res)
  expect_identical(.Random.seed, old_seed)

  set.seed(840840)
  s_best_a <- summary(best_res)
  set.seed(840840)
  s_best_b <- summary(best_res)
  expect_identical(s_best_a, s_best_b)

  # print.summary.best_opt_result should not touch rng and be reproducible
  set.seed(850850)
  old_seed <- .Random.seed
  out_sum2 <- utils::capture.output(print(s_best))
  expect_identical(.Random.seed, old_seed)

  set.seed(860860)
  out_sum_a <- utils::capture.output(print(s_best))
  set.seed(860860)
  out_sum_b <- utils::capture.output(print(s_best))
  expect_identical(out_sum_a, out_sum_b)

  # print.best_opt_result should not touch rng and be reproducible
  set.seed(870870)
  old_seed <- .Random.seed
  out_best2 <- utils::capture.output(print(best_res), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(880880)
  out_best_a <- utils::capture.output(print(best_res), type = "message")
  set.seed(880880)
  out_best_b <- utils::capture.output(print(best_res), type = "message")
  expect_identical(out_best_a, out_best_b)

  # str.best_opt_result should not touch rng and be reproducible
  set.seed(890890)
  old_seed <- .Random.seed
  out_str2 <- utils::capture.output(str(best_res))
  expect_identical(.Random.seed, old_seed)

  set.seed(900900)
  out_str_a <- utils::capture.output(str(best_res))
  set.seed(900900)
  out_str_b <- utils::capture.output(str(best_res))
  expect_identical(out_str_a, out_str_b)

  # plot.best_opt_result should not touch rng and return same object
  set.seed(910910)
  old_seed <- .Random.seed
  res_plot2 <- plot(best_res)
  expect_identical(.Random.seed, old_seed)

  set.seed(920920)
  p1 <- plot(best_res)
  set.seed(920920)
  p2 <- plot(best_res)
  expect_identical(p1, p2)

  # select_opt should not touch rng and be reproducible
  set.seed(930930)
  old_seed <- .Random.seed
  sel_tmp <- select_opt(
    x             = best_res,
    smd_type      = "max",
    smd_groups    = NULL,
    smd_variables = NULL,
    perc_matched  = NULL
  )
  expect_identical(.Random.seed, old_seed)

  set.seed(940940)
  sel_a <- select_opt(
    x             = best_res,
    smd_type      = "max",
    smd_groups    = NULL,
    smd_variables = NULL,
    perc_matched  = NULL
  )
  set.seed(940940)
  sel_b <- select_opt(
    x             = best_res,
    smd_type      = "max",
    smd_groups    = NULL,
    smd_variables = NULL,
    perc_matched  = NULL
  )
  expect_identical(sel_a, sel_b)

  # summary.select_result should not touch rng and be reproducible
  set.seed(950950)
  old_seed <- .Random.seed
  s_sel2 <- summary(sel_res)
  expect_identical(.Random.seed, old_seed)

  set.seed(960960)
  s_sel_a <- summary(sel_res)
  set.seed(960960)
  s_sel_b <- summary(sel_res)
  expect_identical(s_sel_a, s_sel_b)

  # print.summary.select_result should not touch rng and be reproducible
  set.seed(970970)
  old_seed <- .Random.seed
  out_s_sel2 <- utils::capture.output(print(s_sel))
  expect_identical(.Random.seed, old_seed)

  set.seed(980980)
  out_s_sel_a <- utils::capture.output(print(s_sel))
  set.seed(980980)
  out_s_sel_b <- utils::capture.output(print(s_sel))
  expect_identical(out_s_sel_a, out_s_sel_b)

  # print.select_result should not touch rng and be reproducible
  set.seed(990990)
  old_seed <- .Random.seed
  out_sel2 <- utils::capture.output(print(sel_res), type = "message")
  expect_identical(.Random.seed, old_seed)

  set.seed(100010)
  out_sel_a <- utils::capture.output(print(sel_res), type = "message")
  set.seed(100010)
  out_sel_b <- utils::capture.output(print(sel_res), type = "message")
  expect_identical(out_sel_a, out_sel_b)

  # str.select_result should not touch rng and be reproducible
  set.seed(101010)
  old_seed <- .Random.seed
  out_str_sel2 <- utils::capture.output(str(sel_res))
  expect_identical(.Random.seed, old_seed)

  set.seed(102020)
  out_str_sel_a <- utils::capture.output(str(sel_res))
  set.seed(102020)
  out_str_sel_b <- utils::capture.output(str(sel_res))
  expect_identical(out_str_sel_a, out_str_sel_b)

  # get_select_params should not touch rng and be reproducible
  set.seed(103030)
  old_seed <- .Random.seed
  params_all2 <- get_select_params(sel_res)
  expect_identical(.Random.seed, old_seed)

  set.seed(104040)
  pa1 <- get_select_params(sel_res, smd_group = some_group)
  set.seed(104040)
  pa2 <- get_select_params(sel_res, smd_group = some_group)
  expect_identical(pa1, pa2)

  # run_selected_matching should not touch rng and be reproducible
  set.seed(105050)
  old_seed <- .Random.seed
  tmp_match <- run_selected_matching(
    x         = sel_res,
    data      = data,
    formula   = formula,
    smd_group = some_group,
    row       = 1L
  )
  expect_identical(.Random.seed, old_seed)

  set.seed(106060)
  rm1 <- run_selected_matching(
    x         = sel_res,
    data      = data,
    formula   = formula,
    smd_group = some_group,
    row       = 1L
  )
  set.seed(106060)
  rm2 <- run_selected_matching(
    x         = sel_res,
    data      = data,
    formula   = formula,
    smd_group = some_group,
    row       = 1L
  )
  expect_identical(rm1, rm2)
})
