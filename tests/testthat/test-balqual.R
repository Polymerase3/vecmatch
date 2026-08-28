test_that("balqual: argument checks and basic run", {
  # estimate the gps
  gps_matrix <- estimate_gps(status ~ age,
    cancer,
    method = "multinom",
    refernce = "control"
  )

  # drop observations outside the csr
  invisible(capture.output(
    {
      csmatrix <- csregion(gps_matrix)
    },
    file = NULL
  ))

  ## testing a clear run
  withr::with_options(list(warn = -1), {
    # matching the csmatrix
    matched_cancer <- match_gps(csmatrix,
      reference = "control",
      caliper = 1,
      kmeans_cluster = 5
    )
  })

  # basic test run
  expect_no_error(
    invisible(
      balqual(
        matched_cancer,
        status ~ age
      )
    )
  )

  # test mean
  expect_no_error(
    invisible(
      capture.output(
        balqual(matched_cancer, status ~ age, statistic = "mean"),
        file = NULL
      )
    )
  )

  # test max
  expect_no_error(
    invisible(
      capture.output(
        balqual(matched_cancer, status ~ age * sex, statistic = "max"),
        file = NULL
      )
    )
  )

  # break cutoffs
  expect_error(
    invisible(
      capture.output(
        balqual(matched_cancer, status ~ age, cutoffs = 1),
        file = NULL
      )
    ),
    regexp = "length"
  )

  # run balqual once to obtain a quality object
  quality_obj <- balqual(
    matched_cancer,
    status ~ age,
    statistic = "max"
  )

  # capture printed output from str()
  out <- utils::capture.output({
    res <- str(quality_obj)
  })

  # check header string is present somewhere in the output
  expect_true(any(grepl("quality object: matching diagnostics", out,
    fixed = TRUE
  )))

  # str() should return the object invisibly
  expect_identical(res, quality_obj)

  # basic structural sanity checks
  expect_s3_class(quality_obj, "quality")
  expect_true("type" %in% names(quality_obj))
  expect_true("statistic" %in% names(quality_obj))

  # attributes used by str.quality should exist (at least some of them)
  expect_true(!is.null(attr(quality_obj, "original_data_before")) ||
    !is.null(attr(quality_obj, "original_data_after")))

  ## reproducibility tests ----------------------------------------------------
  # 1) balqual() should not change the global .Random.seed
  set.seed(12345)
  old_seed <- .Random.seed

  tmp_quality <- balqual(
    matched_cancer,
    status ~ age,
    statistic = "max"
  )

  expect_identical(.Random.seed, old_seed)

  # 2) running balqual() twice with the same seed should give identical results
  set.seed(4242)
  q1 <- balqual(
    matched_cancer,
    status ~ age,
    statistic = "max"
  )

  set.seed(4242)
  q2 <- balqual(
    matched_cancer,
    status ~ age,
    statistic = "max"
  )

  expect_identical(q1, q2)
})

test_that("balqual: descriptive statistics via type = 'desc_full'", {
  gps_matrix <- estimate_gps(status ~ age * sex,
    cancer,
    method = "multinom",
    reference = "control"
  )

  invisible(capture.output(
    {
      csmatrix <- csregion(gps_matrix)
    },
    file = NULL
  ))

  withr::with_options(list(warn = -1), {
    matched_cancer <- match_gps(csmatrix,
      reference = "control",
      caliper = 1,
      kmeans_cluster = 5
    )
  })

  ## `desc_full` combined with a balance metric -------------------------------
  quality_obj <- balqual(
    matched_cancer,
    status ~ age * sex,
    type = c("smd", "desc_full"),
    statistic = "max",
    cutoffs = 0.2
  )

  desc <- quality_obj$quality_desc

  expect_s3_class(desc, "data.frame")
  expect_identical(
    colnames(desc),
    c(
      "Variable", "Group", "Type", "Statistic", "Before", "After",
      "Percent_Before", "Percent_After"
    )
  )

  # the covariates are described as given in the `formula`: `sex` is described
  # once by its levels rather than once per dummy, and the interaction term is
  # not described at all
  expect_identical(sort(unique(desc$Variable)), c("age", "sex"))
  expect_false(any(grepl(":", desc$Variable, fixed = TRUE)))

  ## the continuous block -----------------------------------------------------
  continuous <- desc[desc$Type == "continuous", ]

  expect_identical(unique(continuous$Variable), "age")

  # the statistics are rows, listed in their canonical order for every group
  n_groups <- nunique(cancer$status)
  expect_identical(
    unique(continuous$Statistic),
    c(
      "N", "Mean", "SD", "Min", "Q1", "Median", "Q3", "Max",
      "Skewness", "Kurtosis"
    )
  )
  expect_identical(nrow(continuous), n_groups * 10L)
  expect_true(all(is.na(continuous$Percent_Before)))
  expect_true(all(is.na(continuous$Percent_After)))
  expect_false(anyNA(continuous$Before))

  # `N` is a statistic like any other, and matches the count table
  count_tab <- as.data.frame(quality_obj$count_table)
  n_row <- continuous[continuous$Statistic == "N" &
    continuous$Group == count_tab$Treatment[1], ]
  expect_identical(as.integer(n_row$Before), as.integer(count_tab$Before[1]))
  expect_identical(as.integer(n_row$After), as.integer(count_tab$After[1]))

  ## the categorical block ----------------------------------------------------
  categorical <- desc[desc$Type == "categorical", ]

  expect_identical(unique(categorical$Variable), "sex")
  expect_identical(unique(categorical$Statistic), levels(cancer$sex))

  # the counts of every level add up to the size of the treatment group ...
  for (stage in c("Before", "After")) {
    totals <- tapply(categorical[[stage]], categorical$Group, sum)
    expect_identical(
      as.integer(totals[count_tab$Treatment[1]]),
      as.integer(count_tab[[stage]][1])
    )
  }

  # ... and the percentages add up to 100 within each group, computed
  # separately for each matching stage
  for (stage in c("Percent_Before", "Percent_After")) {
    perc <- tapply(categorical[[stage]], categorical$Group, sum)
    expect_true(all(abs(perc - 100) < 1e-6))
  }

  # the balance metrics are still reported next to the descriptives
  expect_true(all(quality_obj$quality_max$Coefficient == "SMD"))
  expect_false(is.null(quality_obj$summary_head))

  ## `desc_full` on its own ---------------------------------------------------
  quality_desc_only <- balqual(
    matched_cancer,
    status ~ age * sex,
    type = "desc_full"
  )

  expect_identical(nrow(quality_desc_only$quality_mean), 0L)
  expect_identical(nrow(quality_desc_only$quality_max), 0L)
  expect_false(is.null(quality_desc_only$quality_desc))

  # `statistic` does not apply to the descriptive metrics and is flagged
  expect_warning(
    balqual(
      matched_cancer,
      status ~ age * sex,
      type = "desc_full",
      statistic = "max"
    ),
    regexp = "ignored"
  )

  # but it is not flagged when a balance metric is requested as well
  expect_no_warning(
    balqual(
      matched_cancer,
      status ~ age * sex,
      type = c("smd", "desc_full"),
      statistic = "max",
      cutoffs = 0.2
    )
  )

  ## the descriptive metrics are not counted when validating the cutoffs ------
  expect_no_error(
    balqual(matched_cancer, status ~ age,
      type = c("smd", "desc_full"),
      cutoffs = 0.2
    )
  )

  expect_error(
    balqual(matched_cancer, status ~ age,
      type = c("smd", "desc_full"),
      cutoffs = c(0.1, 0.2)
    ),
    regexp = "length"
  )

  ## the default call is unaffected -------------------------------------------
  expect_null(balqual(matched_cancer, status ~ age)$quality_desc)

  ## `desc_reduced` is the same table restricted to four statistics -----------
  quality_reduced <- balqual(
    matched_cancer,
    status ~ age * sex,
    type = c("smd", "desc_reduced"),
    statistic = "max",
    cutoffs = 0.2
  )

  reduced <- quality_reduced$quality_desc

  expect_identical(colnames(reduced), colnames(desc))

  # only the continuous statistics are reduced ...
  expect_identical(
    unique(reduced$Statistic[reduced$Type == "continuous"]),
    c("N", "Min", "Mean", "Median", "Max")
  )

  # ... and the rows that remain hold the same values as in `desc_full`. The
  # two differ only in the order of the statistics, so both are sorted first.
  normalize <- function(df) {
    df <- df[order(df$Variable, df$Group, df$Statistic), ]
    rownames(df) <- NULL
    df
  }

  expect_identical(
    normalize(reduced),
    normalize(desc[
      desc$Type == "categorical" |
        desc$Statistic %in% c("N", "Min", "Mean", "Median", "Max"),
    ])
  )

  # the cross-tabulation of the categorical covariates is untouched
  expect_identical(
    normalize(reduced[reduced$Type == "categorical", ]),
    normalize(desc[desc$Type == "categorical", ])
  )

  # the two descriptive metrics are mutually exclusive
  expect_error(
    balqual(matched_cancer, status ~ age,
      type = c("desc_full", "desc_reduced")
    ),
    regexp = "mutually exclusive"
  )

  expect_warning(
    balqual(matched_cancer, status ~ age * sex,
      type = "desc_reduced",
      statistic = "mean"
    ),
    regexp = "ignored"
  )

  ## printing and str() include the descriptive table -------------------------
  out <- utils::capture.output(print(quality_obj))
  expect_true(any(grepl("Descriptive statistics of the continuous covariates",
    out,
    fixed = TRUE
  )))
  expect_true(any(grepl("Distribution of the categorical covariates", out,
    fixed = TRUE
  )))

  # the counts and percentages are compressed into a single `N (%)` cell
  expect_true(any(grepl("175 (46.9%)", out, fixed = TRUE)))

  # a formula without categorical covariates prints no categorical block
  out_numeric <- utils::capture.output(
    print(balqual(matched_cancer, status ~ age, type = "desc_full"))
  )
  expect_true(any(grepl("continuous covariates", out_numeric, fixed = TRUE)))
  expect_false(any(grepl("categorical covariates", out_numeric, fixed = TRUE)))

  out_str <- utils::capture.output(str(quality_obj))
  expect_true(any(grepl("quality_desc:", out_str, fixed = TRUE)))
})

test_that("balqual: cutoffs are matched to their metric by name", {
  gps_matrix <- estimate_gps(status ~ age,
    cancer,
    method = "multinom",
    reference = "control"
  )

  invisible(capture.output(
    {
      csmatrix <- csregion(gps_matrix)
    },
    file = NULL
  ))

  withr::with_options(list(warn = -1), {
    matched_cancer <- match_gps(csmatrix,
      reference = "control",
      caliper = 1,
      kmeans_cluster = 5
    )
  })

  # `var_ratio` must be judged against its own cutoff (2), not against the
  # `smd` cutoff that precedes it in the `type` argument
  expect_no_warning(
    quality_obj <- balqual(
      matched_cancer,
      status ~ age,
      type = c("smd", "var_ratio"),
      statistic = "max",
      cutoffs = c(0.1, 2)
    )
  )

  var_rows <- quality_obj$quality_max[
    quality_obj$quality_max$Coefficient == "Var",
  ]

  expect_identical(
    var_rows$Reduction,
    ifelse(var_rows$After < 2, "Balanced", "Not Balanced")
  )
})

test_that(".desc_stats_table: moment estimators and NA handling", {
  x <- c(1, 2, 3, 4, 10)
  covs <- matrix(x, ncol = 1, dimnames = list(NULL, "v"))

  res <- .desc_stats_table(covs, rep("g", 5), time = "Before", round = 5)

  expect_identical(res$N, 5L)
  expect_equal(res$Mean, 4)
  expect_equal(res$Median, 3)
  expect_equal(res$Min, 1)
  expect_equal(res$Max, 10)
  expect_equal(res$Q1, 2)
  expect_equal(res$Q3, 4)

  # classical (type 1) moment estimators
  expect_equal(res$Skewness, round(36 / 10^1.5, 5))
  expect_equal(res$Kurtosis, round(278.8 / 100 - 3, 5))

  # constant covariates have no defined shape statistics
  const <- matrix(rep(1, 5), ncol = 1, dimnames = list(NULL, "v"))
  res_const <- .desc_stats_table(const, rep("g", 5), time = "Before")

  expect_identical(res_const$SD, 0)
  expect_true(is.na(res_const$Skewness))
  expect_true(is.na(res_const$Kurtosis))

  # missing values are dropped per covariate
  with_na <- matrix(c(1, 2, NA, 4, 5), ncol = 1, dimnames = list(NULL, "v"))
  res_na <- .desc_stats_table(with_na, rep("g", 5), time = "After")

  expect_identical(res_na$N, 4L)
  expect_equal(res_na$Mean, 3)
})

test_that(".desc_table: dispatches covariates by their type", {
  covs <- data.frame(
    num = c(1, 2, 3, 4, 10, 20),
    fct = factor(c("b", "a", "a", "b", "a", "a"), levels = c("b", "a")),
    chr = c("x", "x", "y", "y", "y", "y"),
    lgl = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  treat <- factor(rep(c("t1", "t0"), each = 3), levels = c("t1", "t0"))

  res <- .desc_table(covs, treat, time = "Before", round = 4)

  # numeric covariates are summarized, everything else cross-tabulated
  expect_identical(
    unique(res$Type[res$Variable == "num"]), "continuous"
  )
  expect_identical(
    sort(unique(res$Variable[res$Type == "categorical"])),
    c("chr", "fct", "lgl")
  )

  # the level order of a factor is preserved, characters are sorted
  expect_identical(unique(res$Level[res$Variable == "fct"]), c("b", "a"))
  expect_identical(unique(res$Level[res$Variable == "chr"]), c("x", "y"))

  # counts and percentages are computed within a treatment level
  fct_t1 <- res[res$Variable == "fct" & res$Group == "t1", ]
  expect_identical(fct_t1$N, c(1L, 2L))
  expect_equal(fct_t1$Percent, c(100 / 3, 200 / 3), tolerance = 1e-3)

  # the treatment level order of a factor is kept
  expect_identical(unique(res$Group), c("t1", "t0"))

  # missing values are excluded from both the counts and the denominator
  covs_na <- data.frame(fct = factor(c("a", "b", NA, "a")))
  res_na <- .desc_table(covs_na, rep("g", 4), time = "After")

  expect_identical(res_na$N, c(2L, 1L))
  expect_equal(res_na$Percent, c(200 / 3, 100 / 3), tolerance = 1e-3)
})

test_that(".desc_to_wide: places the matching stages side by side", {
  covs <- data.frame(
    num = c(1, 2, 3, 4, 10, 20),
    fct = factor(c("a", "b", "c", "a", "b", "c")),
    stringsAsFactors = FALSE
  )
  treat <- rep("g", 6)

  # after matching, level "c" no longer occurs
  kept <- covs$fct != "c"

  desc_long <- rbind(
    .desc_table(covs, treat, time = "Before", round = 4),
    .desc_table(covs[kept, ], treat[kept], time = "After", round = 4)
  )

  out <- .desc_to_wide(
    desc_long,
    stats_keep = c("N", "Mean", "Max"),
    var_order = names(covs),
    group_order = "g"
  )

  expect_identical(
    colnames(out),
    c(
      "Variable", "Group", "Type", "Statistic", "Before", "After",
      "Percent_Before", "Percent_After"
    )
  )

  # the requested statistics become rows, in the order they were requested
  continuous <- out[out$Type == "continuous", ]
  expect_identical(continuous$Statistic, c("N", "Mean", "Max"))
  expect_identical(continuous$Before, c(6, round(40 / 6, 4), 20))
  expect_true(all(is.na(continuous$Percent_Before)))

  # a level dropped entirely by matching is reported as observed zero times,
  # not as missing
  dropped <- out[out$Type == "categorical" & out$Statistic == "c", ]
  expect_identical(dropped$Before, 2)
  expect_identical(dropped$After, 0)
  expect_identical(dropped$Percent_After, 0)

  # the percentages are computed within each matching stage separately, and
  # add up to 100 up to the rounding requested through `round`
  categorical <- out[out$Type == "categorical", ]
  expect_equal(sum(categorical$Percent_Before), 100, tolerance = 1e-3)
  expect_equal(sum(categorical$Percent_After), 100, tolerance = 1e-3)
})
