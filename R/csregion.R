#' @title Filter the Data Based on Common Support Region
#'
#' @description The `csregion()` function estimates the boundaries of the
#'   rectangular common support region, as defined by Lopez and Gutman (2017),
#'   and filters the matrix of generalized propensity scores based on these
#'   boundaries. The function returns a matrix of observations whose generalized
#'   propensity scores lie within the treatment group-specific boundaries.
#'
#' @param gps_matrix An object of classes `gps` and `data.frame` (e.g., created
#'   by the `estimate_gps()` function). The first column corresponds to the
#'   treatment or grouping variable, while the other columns represent the
#'   treatment assignment probabilities calculated separately for each
#'   hypotetical treatment group. The number of columns should therefore be
#'   equal to the number of unique levels of the treatment variable plus one
#'   (for the treatment variable itself). The number of rows should correspond
#'   to the number of subjects for which generalized propensity scores were
#'   estimated.
#'
#' @param borders A character string specifying how to handle observations at
#'   the edges of the Common Support Region (CSR). Acceptable values are
#'   `"include"` and `"exclude"`. If `"include"` is selected (default),
#'   observations with Generalized Propensity Scores (GPS) exactly equal to the
#'   CSR boundaries are retained for further analysis. This corresponds to a
#'   non-strict inequality: \code{lower_bound <= GPS <= upper_bound}. If
#'   `"exclude"` is selected, observations lying exactly on the CSR boundaries
#'   are removed. This corresponds to a strict inequality: \code{lower_bound <
#'   GPS < upper_bound}. Using `"exclude"` will typically result in a slightly
#'   smaller matched sample size compared to `"include"`, but may be preferred
#'   for more conservative matching.
#'
#' @param refit Logical. If `TRUE` (default), the model used to estimate the GPS
#'   is refitted after excluding samples outside the common support region,
#'   using the same formula and method as in the original `estimate_gps()` call.
#'   If `FALSE`, the model is not refitted, but still only samples within the
#'   CSR are retained. Refitting is recommended, as suggested by Lopez and
#'   Gutman (2017).
#'
#' @return A numeric matrix similar to the one returned by `estimate_gps()`,
#'   but with the number of rows reduced to exclude those observations that do
#'   not fit within the common support region (CSR) boundaries. The returned
#'   object also possesses additional attributes that summarize the calculation
#'   process of the CSR boundaries:
#'  * `filter_matrix` - A logical matrix with the same dimensions as the
#'  gps-part of `gps_matrix`, indicating which treatment assignment
#'  probabilities fall within the CSR boundaries,
#'  * `filter_vector` - A vector indicating whether each observation was kept
#'  (`TRUE`) or removed (`FALSE`), essentially a row-wise
#'  sum of `filter_matrix`,
#'  * `csr_summary` - A summary of the CSR calculation process, including
#'  details of the boundaries and the number of observations filtered.
#'  * `csr_data` - The original dataset used for the estimation of generalized
#'  propensity scores (`original_data` attribute of the `gps` object) filtered
#'  by the `filter_vector`
#'
#' @examples
#' # We could estimate simples generalized propensity scores for the `iris`
#' # dataset
#' gps <- estimate_gps(Species ~ Sepal.Length, data = iris)
#'
#' # And then define the common support region boundaries using `csregion()`
#' gps_csr <- csregion(gps)
#'
#' # The additional information of the CSR-calculation process are
#' # accessible through the attributes described in the `*Value*` section
#' attr(gps_csr, "filter_matrix")
#' attr(gps_csr, "csr_summary")
#' attr(gps_csr, "csr_data")
#'
#' @export
csregion <- function(gps_matrix,
                     borders = "include",
                     refit = TRUE) {
  csr_data <- attr(gps_matrix, "original_data")

  .chk_cond(
    "gps" %nin% class(gps_matrix),
    "The `gps_matrix` argument must be of class `gps`."
  )

  # check borders arg
  chk::chk_character(borders)
  chk::chk_length(borders, length = 1)
  .chk_cond(
    borders %nin% c("include", "exclude"),
    'The `borders` argument can only take one of the
            following values: "include", "exclude".'
  )

  # check refit arg
  chk::chk_flag(refit)

  ## Calculating the csr_low
  csr_low <- apply(
    stats::aggregate(. ~ treatment,
      data = gps_matrix,
      FUN = function(x) min(x)
    )[, -1],
    2,
    max
  )

  ## Calculating the csr_high
  csr_high <- apply(
    stats::aggregate(. ~ treatment,
      data = gps_matrix,
      FUN = function(x) max(x)
    )[, -1],
    2,
    min
  )

  ## filter out the unvalid observations
  filter_matrix <- mapply(function(df_col, vec_low, vec_high) {
    switch(borders,
      include = vec_low <= df_col & df_col <= vec_high,
      exclude = vec_low < df_col & df_col < vec_high
    )
  }, gps_matrix[, 2:ncol(gps_matrix)], csr_low, csr_high)

  ## summarizing the logical matrix to a subset vector
  filter_vector <- apply(filter_matrix, 1, all)

  ## defining the number of negatives
  n_negative <- sum(!filter_vector)
  n_negative_matrix <- colSums(!filter_matrix)

  ## refitting the gps_matrix
  if (refit) {
    # save original gps_matrix restricted to CSR as a fallback
    gps_matrix_csr <- subset(gps_matrix, filter_vector)

    # limit original data only to the CSR
    csr_filtered <- csr_data[filter_vector, , drop = FALSE]

    # prepare function call for estimate_gps
    estimate_call <- attr(gps_matrix, "function_call")
    estimate_call$data <- csr_filtered

    # try to evaluate the changed call
    gps_matrix_refit <- tryCatch(
      eval(estimate_call),
      error = function(e) {
        chk::wrn(strwrap(
          "Refitting of the GPS model on the CSR-restricted data failed. ",
          "Falling back to the original GPS restricted to the common support ",
          "region (refit = FALSE).",
          prefix = " ",
          initial = ""
        ))
        NULL
      }
    )

    if (is.null(gps_matrix_refit)) {
      # refitting failed: use the non-refitted CSR-restricted gps_matrix
      gps_matrix <- gps_matrix_csr
      refit <- FALSE # only local, used if you want to record this below
    } else {
      # refitting succeeded
      gps_matrix <- gps_matrix_refit
      refit <- TRUE
    }
  } else {
    ## subset the gps_matrix without refitting
    gps_matrix <- subset(gps_matrix, filter_vector)
  }

  ## drop unused levels of treatment variable
  gps_matrix[, "treatment"] <- droplevels(gps_matrix[, "treatment"])

  ## detect low number of observations in groups and print a warning
  if (any(table(gps_matrix[, "treatment"]) < 20)) {
    chk::wrn(strwrap("Some groups have fewer than 20 observations, which may
    impact the performance of the matching process. Consider using
    `replace = TRUE`in `match_gps()` to address this.",
      prefix = " ", initial = ""
    ))
  }

  ## assembling the print results list
  csr_summary <- data.frame(
    treatment = colnames(gps_matrix)[2:ncol(gps_matrix)],
    csr_low = csr_low,
    csr_high = csr_high,
    n_negative_matrix = n_negative_matrix
  )

  ## add attributes and set class
  gps_matrix <- structure(
    gps_matrix,
    filter_matrix = filter_matrix,
    filter_vector = filter_vector,
    csr_summary   = csr_summary,
    csr_data      = csr_data[filter_vector, ],
    class         = c("csr", "gps", "data.frame")
  )

  # return the gps_matrix
  return(gps_matrix)
}

#' @export
summary.csr <- function(object, digits = 3, ...) {
  csr_summary <- attr(object, "csr_summary")
  filter_vector <- attr(object, "filter_vector")

  if (is.null(csr_summary) || is.null(filter_vector)) {
    cli::cli_abort(
      "Object of class {.cls csr} is missing required attributes
       {.field csr_summary} or {.field filter_vector}."
    )
  }

  n_total <- length(filter_vector)
  n_kept <- sum(filter_vector)
  n_excluded <- n_total - n_kept

  # sizes AFTER filtering (in actual csr object)
  group_sizes_after <- table(object[["treatment"]])

  # align to treatment order in csr_summary
  tr <- csr_summary$treatment
  after <- as.integer(group_sizes_after[tr])
  after[is.na(after)] <- 0L
  before <- after + csr_summary$n_negative_matrix

  # store numeric tables, formatting happens in print.summary.csr()
  # enables programmatic use of summary.csr()
  df_bounds <- data.frame(
    Treatment = tr,
    csr_low = csr_summary$csr_low,
    csr_high = csr_summary$csr_high,
    n_excluded = csr_summary$n_negative_matrix,
    check.names = FALSE
  )

  df_groups <- data.frame(
    Treatment = tr,
    n_before = before,
    n_after = after,
    Excluded = csr_summary$n_negative_matrix,
    check.names = FALSE
  )

  res <- list(
    csr_summary        = csr_summary,
    filter_vector      = filter_vector,
    n_total            = n_total,
    n_kept             = n_kept,
    n_excluded         = n_excluded,
    group_sizes_before = stats::setNames(before, tr),
    group_sizes_after  = stats::setNames(after, tr),
    bounds_table       = df_bounds,
    groups_table       = df_groups,
    digits             = digits
  )

  class(res) <- "summary.csr"
  return(res)
}

#' @export
print.summary.csr <- function(x, digits = NULL, ...) {
  # decide digits: explicit argument > stored in object > default 3
  if (is.null(digits)) {
    digits <- if (!is.null(x$digits)) x$digits else 3L
  }

  fmt <- function(z) format(round(z, digits = digits), trim = TRUE, nsmall = digits)

  # local helper: fixed-width table; may convert in the future to a CLI-like
  # table
  print_table <- function(df) {
    colnames_df <- colnames(df)

    cat("\n")
    cat(paste(sprintf("%-18s", colnames_df), collapse = " | "), "\n")
    cat(paste(rep("-", 20 * ncol(df)), collapse = ""), "\n")

    apply(df, 1, function(row) {
      cat(paste(sprintf("%-18s", row), collapse = " | "), "\n")
    })
    cat("\n")
  }

  # build printable tables from stored numeric tables
  bounds_raw <- x$bounds_table
  groups_raw <- x$groups_table

  df_bounds <- data.frame(
    Treatment = bounds_raw$Treatment,
    "Lower CSR limit" = fmt(bounds_raw$csr_low),
    "Upper CSR limit" = fmt(bounds_raw$csr_high),
    "Number excluded" = bounds_raw$n_excluded,
    check.names = FALSE
  )

  df_groups <- data.frame(
    Treatment = groups_raw$Treatment,
    "N before" = groups_raw$n_before,
    "N after" = groups_raw$n_after,
    Excluded = groups_raw$Excluded,
    check.names = FALSE
  )

  # ---------------- CLI output ----------------
  cli::cli_h1("Rectangular CSR Borders Evaluation")
  cli::cli_text("Class: {.cls csr} (gps filtered to common support region)")
  cli::cli_text(
    "Rows kept: {x$n_kept} / {x$n_total} (excluded: {x$n_excluded})"
  )
  cli::cli_text(
    "Columns: {length(x$group_sizes_after)} treatment levels (plus GPS columns in the original object)"
  )
  cli::cli_rule()

  cli::cli_h2("Per-treatment CSR bounds")
  print_table(df_bounds)

  cli::cli_h2("Group sizes before and after CSR filtering")
  print_table(df_groups)

  cli::cli_rule()
  cli::cli_text(
    "Details of the CSR calculation can also be accessed via ",
    "{.code attr(x, 'csr_summary')}, ",
    "{.code attr(x, 'filter_matrix')}, ",
    "and {.code attr(x, 'csr_data')} on the original {.cls csr} object."
  )

  invisible(x)
}

#' @export
print.csr <- function(x, ...) {
  cli::cli_text("{.strong csr object} (gps filtered to common support region)")
  .print_gps_core(x, ...) # helper defined in estimate_gps.R
}

#' @export
#' @export
str.csr <- function(object, ...) {
  n <- nrow(object)
  p <- ncol(object)

  treatment <- object[["treatment"]]

  trt_levels <- levels(treatment)

  group_sizes_after <- if (!is.null(treatment)) table(treatment) else NULL

  filter_vector <- attr(object, "filter_vector")
  csr_summary <- attr(object, "csr_summary")
  csr_data <- attr(object, "csr_data")
  filter_matrix <- attr(object, "filter_matrix")

  # totals
  n_total <- length(filter_vector)
  n_kept <- sum(filter_vector)
  n_excluded <- n_total - n_kept


  # ---- header -------------------------------------------------------------
  cat("csr object: gps filtered to common support region\n")
  cat(sprintf(" Dimensions (csr object): %d rows x %d columns\n", n, p))
  cat(sprintf(" Treatment column: %s\n", "treatment"))
  cat(sprintf(
    " Treatment levels: %s\n",
    if (length(trt_levels)) paste(trt_levels, collapse = ", ") else "<none>"
  ))
  if (!is.null(group_sizes_after)) {
    cat(" Group sizes after CSR filtering:\n")
    gs_txt <- paste(
      sprintf("  - %s: %d", names(group_sizes_after), as.integer(group_sizes_after)),
      collapse = "\n"
    )
    cat(gs_txt, "\n", sep = "")
  }

  cat(sprintf(
    " Rows kept: %d / %d (excluded: %d)\n",
    n_kept, n_total, n_excluded
  ))

  cat(sprintf(
    " csr_data: data.frame with %d rows and %d columns\n",
    NROW(csr_data), NCOL(csr_data)
  ))

  cat(sprintf(
    " filter_matrix: logical matrix %d x %d\n",
    NROW(filter_matrix), NCOL(filter_matrix)
  ))

  cat(sprintf(
    " csr_summary: data.frame with %d rows (per-treatment bounds)\n",
    NROW(csr_summary)
  ))

  cat("\nUnderlying data.frame structure:\n")

  # ---- delegate to data.frame str() ---------------------------------------
  utils::str(unclass(object), ...)

  invisible(object)
}

#' @export
as.data.frame.csr <- function(x, ...) {
  class(x) <- "data.frame"
  NextMethod()
}

#' @export
plot.csr <- function(x,
                     gps_cols = NULL,
                     treatment_col = "treatment",
                     ...) {
  # x: csr object (data.frame with GPS columns + treatment)
  # gps_cols: optional character vector of GPS columns to plot (panels)
  # treatment_col: name of treatment column in x
  # ...: additional arguments passed to boxplot()

  # --- basic checks ----------------------------------------------------------
  df <- as.data.frame(x)

  trt <- df[[treatment_col]]
  trt <- factor(trt)

  all_gps_cols <- setdiff(names(df), treatment_col)

  if (is.null(gps_cols)) {
    gps_cols <- all_gps_cols
  } else {
    gps_cols <- intersect(gps_cols, all_gps_cols)
    if (!length(gps_cols)) {
      stop(
        "None of the requested `gps_cols` are valid. ",
        "Available: ", paste(all_gps_cols, collapse = ", ")
      )
    }
  }

  n_gps <- length(gps_cols)

  # --- layout: number of panels based on number of GPS columns --------------
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  if (n_gps == 1L) {
    graphics::par(mfrow = c(1, 1))
  } else if (n_gps == 2L) {
    graphics::par(mfrow = c(1, 2))
  } else if (n_gps <= 4L) {
    graphics::par(mfrow = c(2, 2))
  } else {
    nrow <- ceiling(n_gps / 2)
    graphics::par(mfrow = c(nrow, 2))
  }

  # --- draw boxplots: one panel per GPS column -------------------------------
  for (col in gps_cols) {
    y <- df[[col]]

    graphics::boxplot(
      y ~ trt,
      main = paste("GPS column:", col),
      xlab = "Treatment",
      ylab = "Generalized propensity score",
      ...
    )
  }

  invisible(x)
}
