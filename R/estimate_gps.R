#' @title Calculate Treatment Allocation Probabilities
#'
#' @description `estimate_gps()` computes generalized propensity scores for
#'   treatment groups by applying a user-defined formula and method. It returns
#'   a matrix of GPS probabilities for each subject and treatment group
#'
#' @param formula a valid R formula, which describes the model used to
#'   calculating the probabilities of receiving a treatment. The variable to be
#'   balanced is on the left side, while the covariates used to predict the
#'   treatment variable are on the right side. To define the interactions
#'   between covariates, use `*`. For more details, refer to [stats::formula()].
#' @param data a data frame with columns specified in the `formula` argument.
#' @param method a single string describing the model used for the calculation
#'   of generalized propensity scores. The default value is set to `multinom`.
#'   For available methods refer to the Details section below.
#' @param link a single string; determines an alternative model for a method
#'   used for estimation. For available links, see Details.
#' @param subset a logical atomic vector of length equal to the number of rows
#'   in the `data` arguments. Allows to filter out observations from the further
#'   analysis, for which the value of the vector is equal to `FALSE`.
#' @param reference a single string describing one class from the treatment
#'   variable, referred to as the baseline category in the calculation of
#'   generalized propensity scores.
#' @param by a single string with the name of a column, contained in the `data`
#'   argument. The dataset will be divided by the groups created by the grouping
#'   `by` variable and the calculation of the propensity scores will be carried
#'   out separately for each group. The results will then be merged and
#'   presented to the user as a single GPS matrix.
#' @param ordinal_treat an atomic vector of the length equal to the length of
#'   unique levels of the treatment variable. Confirms, that the treatment
#'   variable is an ordinal variable and adjusts its levels, to the order of
#'   levels specified in the argument. Is a call to the function
#'   `factor(treat, levels = ordinal_treat, ordered = TRUE`.
#' @param fit_object a logical flag. If `TRUE`, the the fitted object is
#'   returned instead of the GPS matrix.
#' @param verbose_output a logical flag. If `TRUE` a more verbose version of the
#'   function is run and the output is printed out to the console.
#' @param ... additional arguments, that can be passed to the fitting function
#'   and are not controlled by the above arguments. For more details and
#'   examples refer to the Details section and documentations of corresponding
#'   functions.
#'
#' @return A data frame of class `gps` with the number of columns equal to
#'   the number of unique treatment variable levels plus one (for the treatment
#'   variable itself) and the number of rows equal to the number of subjects in
#'   the initial dataset. The original dataset used for estimation can be
#'   accessed as the `original_data` attribute.
#'
#' @srrstats {G1.3} Key statistical terms used in `vecmatch`, including
#'   generalized propensity scores (GPS), vector matching, common support
#'   region, and treatment allocation probabilities, are defined and explained
#'   in the documentation of `estimate_gps()`, `csregion()`, `match_gps()`
#' @srrstats {G1.4} roxygen2 is used to document all functions
#'
#' @details The main goal of the `estimate_gps()` function is to calculate the
#'   generalized propensity scores aka. treatment allocation probabilities. It
#'   is the first step in the workflow of the vector matching algorithm and is
#'   essential for the further analysis. The returned matrix of class `gps` can
#'   then be passed to the `csregion()` function to calculate the rectangular
#'   common support region boundaries and drop samples not eligible for the
#'   further analysis. The list of available methods operated by the
#'   `estimate_gps()` is provided below with a short description and function
#'   used for the calculations:
#'   * `multinom` - multinomial logistic regression model [nnet::multinom()]
#'   * `vglm` - vector generalized linear model for multinomial data
#'   [VGAM::vglm()],
#'   * `brglm2` - bias reduction model for multinomial responses using the
#'   Poisson trick [brglm2::brmultinom()],
#'   * `mblogit` - baseline-category logit models [mclogit::mblogit()].
#'   * `polr` - ordered logistic or probit regression only for ordered factor
#'   variables from [MASS::polr()]. The `method` argument of the underlying
#'   `MASS::polr()` package function can be controlled with the `link` argument.
#'   Available options: `link = c("logistic", "probit", "loglog", "cloglog",
#'   "cauchit")`
#'
#' @seealso [csregion()] for the calculation of common support region,
#'   [match_gps()] for the matching of generalized propensity scores
#' @examples
#'
#' ## Example 1: multinomial bias-reduced model (brglm2) on `airquality`
#' if (requireNamespace("brglm2", quietly = TRUE)) {
#'   library(brglm2)
#'
#'   # Initial imbalance of means
#'   tapply(airquality$Wind, airquality$Month, mean, na.rm = TRUE)
#'
#'   # Formula definition
#'   formula_air <- Month ~ Wind
#'
#'   # Estimating the generalized propensity scores using brglm2
#'   gp_scores <- estimate_gps(
#'     formula_air,
#'     data = airquality,
#'     method = "brglm2",
#'     reference = "5",
#'     verbose_output = TRUE,
#'     control = brglm2::brglmControl(type = "MPL_Jeffreys")
#'   )
#'
#'   # Filtering the observations outside the csr region
#'   gps_csr <- csregion(gp_scores)
#'   filter_which <- attr(gps_csr, "filter_vector")
#'   filtered_air <- airquality[filter_which, ]
#'
#'   # Imbalance after csr
#'   tapply(filtered_air$Wind, filtered_air$Month, mean, na.rm = TRUE)
#'
#'   # Visual inspection using raincloud
#'   raincloud(
#'     filtered_air,
#'     y = Wind,
#'     group = Month,
#'     significance = "t_test"
#'   )
#' }
#'
#' ## Example 2: ordered treatment, subset, by, and non-default link
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   library(MASS)
#'
#'   # Prepare a clean subset of `airquality`
#'   aq <- subset(
#'     airquality,
#'     !is.na(Month) & !is.na(Wind) & !is.na(Temp)
#'   )
#'
#'   # Grouping variable used in the `by` argument
#'   aq$Temp_group <- ifelse(
#'     aq$Temp > stats::median(aq$Temp, na.rm = TRUE),
#'     "high",
#'     "low"
#'   )
#'
#'   # Explicit order for the (ordinal) treatment variable
#'   ord_month <- sort(unique(aq$Month))
#'
#'   gps_ord <- estimate_gps(
#'     Month ~ Wind + Temp,
#'     data = aq,
#'     method = "polr",
#'     link = "probit",
#'     subset = aq$Wind > stats::median(aq$Wind, na.rm = TRUE),
#'     by = "Temp_group",
#'     ordinal_treat = ord_month,
#'     reference = "5"
#'   )
#' }
#'
#' @export
estimate_gps <- function(formula,
                         data = NULL,
                         method = "multinom",
                         link = NULL,
                         reference = NULL,
                         by = NULL,
                         subset = NULL,
                         ordinal_treat = NULL,
                         fit_object = FALSE,
                         verbose_output = FALSE,
                         ...) {
  ####################### INPUT CHECKING  ######################################
  #######################       AND       ######################################
  ####################### DATA PROCESSING ######################################
  call <- match.call()
  orig_formula <- formula # store the original as well (eg. when formula is an
  #                         object)
  call[["formula"]] <- orig_formula
  args <- list(...)

  dots <- substitute(list(...))[-1]

  ## If function call then substitute, else normal
  additional_args <- which(names(call)[-1] %nin% names(formals(estimate_gps)))

  if (!is.null(dots)) {
    for (i in seq_along(additional_args)) {
      argname <- names(args)[i]
      callname <- paste0(argname, "_call")
      callname_char <- paste0(callname, "_char")
      args[callname] <- FALSE

      if (is.call(dots[[i]])) {
        args[callname] <- TRUE
        args[callname_char] <- as.character(dots[[i]])[1]
      }
    }
  }

  # formula
  data_list <- .process_formula(formula, data)

  # args assignment to list used in calculations
  args["treat"] <- list(data_list[["treat"]])
  args["covs"] <- list(data_list[["model_covs"]])

  # process and check ordinal_treat
  if (!is.null(ordinal_treat)) {
    chk::chk_atomic(ordinal_treat)
    chk::chk_vector(ordinal_treat)

    .chk_cond(
      length(ordinal_treat) != length(unique(args[["treat"]])),
      "The numbers of levels provided in `ordinal_treat` has to
                     be the same, as the number of unique levels in the
                     treatment variable."
    )

    args[["treat"]] <- factor(args[["treat"]],
      levels = ordinal_treat,
      ordered = TRUE
    )
  } else {
    args[["treat"]] <- factor(args[["treat"]],
      levels = unique(args[["treat"]]),
      ordered = FALSE
    )
  }

  # data
  if (!is.null(data)) .check_df(data)

  # method
  if (is.null(substitute(method)) || missing(method)) {
    method <- "multinom"
    attr(method, "name") <- method
  } else {
    .check_method(method)
    method_name <- deparse1(substitute(method))
    attr(method, "name") <- method_name
  }

  # link
  available_links <- .gps_methods[[method]]$link_fun

  if (!is.null(link)) {
    chk::chk_string(link)

    .chk_cond(
      link %nin% available_links,
      sprintf(
        "The argument `link` for the method %s only accepts values: %s",
        method, word_list(add_quotes(available_links))
      )
    )

    args[["link"]] <- link
  } else {
    args[["link"]] <- available_links[1]
  }

  # reference
  ref_list <- .process_ref(args[["treat"]],
    ordinal_treat = ordinal_treat,
    reference = reference
  )

  args[["treat"]] <- ref_list[["data.relevel"]]
  reference <- ref_list[["reference"]]

  # subset
  if (!is.null(subset)) {
    chk::chk_string(subset)

    .chk_cond(
      subset %nin% colnames(data),
      sprintf(
        "The column %s defined in the `subset` argument was not found in
                the provided dataset.", subset
      )
    )


    subset_logvec <- as.vector(data[[subset]])

    .chk_cond(
      !is.logical(subset_logvec) || length(dim(subset_logvec)) == 2L,
      "The `subset` argument has to be a name of single column with
              logical values."
    )

    use_subset <- TRUE
  } else {
    use_subset <- FALSE
  }

  # fit_object and verbose_output check
  chk::chk_all(list(fit_object, verbose_output), chk::chk_flag)

  # assembling the arguments list
  if (use_subset) {
    args[["treat"]] <- args[["treat"]][subset_logvec]
    args["covs"] <- list(
      data_list[["reported_covs"]][subset_logvec, , drop = FALSE]
    )

    args[".data"] <- list(data[subset_logvec, , drop = FALSE])
    args[["by"]] <- .process_by(by, data, args[["treat"]])[subset_logvec, ,
      drop = FALSE
    ]
  } else {
    args["covs"] <- list(data_list[["reported_covs"]])
    args[".data"] <- list(data)
    args[["by"]] <- .process_by(by, data, args[["treat"]])
  }

  args["formula"] <- list(formula)
  args["method"] <- list(method)
  args["reference"] <- reference
  args["fit_object"] <- list(fit_object)
  args["verbose_output"] <- list(verbose_output)

  ####################### FITTING ##############################################
  fit_func <- .gps_methods[[method]]$func_used
  use_by <- FALSE
  if (!is.null(args[["by"]])) {
    fitted_object <- list()
    treat_by <- list()
    by.levels <- levels(attr(args[["by"]], "by.factor"))
    use_by <- TRUE

    for (i in seq_along(by.levels)) {
      # subset rule
      by_sub <- attr(args[["by"]], "by.factor") == by.levels[i]

      # create env and subset vars
      by_env <- list2env(args, envir = new.env(), parent = emptyenv())

      with(by_env, {
        selected <- mget(c(".data", "covs", "treat"), envir = by_env)
        subsetted <- lapply(selected, function(x) {
          if (is.data.frame(x)) {
            x[by_sub, , drop = FALSE]
          } else if (is.atomic(x)) {
            x[by_sub]
          }
        })
      })

      # overwrite
      list2env(by_env$subsetted, envir = by_env)

      withr::with_preserve_seed({
        # model the data
        fit <- do.call(
          fit_func,
          as.list(by_env)
        )
      })

      # append to existing list
      fitted_object <- append(fitted_object, list(fit))

      # delete env
      rm(by_env)

      ## save treatment var
      treat_by[[i]] <- args[["treat"]][by_sub]
    }
  } else {
    withr::with_preserve_seed({
      fitted_object <- do.call(
        fit_func,
        args
      )
    })
  }

  ####################### OUTPUT OBJECT ########################################
  # Define the class of the output
  if (fit_object) {
    return(fitted_object)
  }

  results <- if (use_by) {
    # Handle case where `use_by` is TRUE
    fitted_values <- if (all(vapply(
      fitted_object, isS4,
      logical(1L)
    )) && method == "vglm") {
      lapply(fitted_object, VGAM::fitted.values)
    } else {
      lapply(fitted_object, "[[", "fitted.values")
    }

    # saving the primary levels, to reset it later
    # cbind converts factors to numeric by default
    treat_labels <- levels(unlist(treat_by))
    treat_levels <- sort(unique(as.integer(unlist(treat_by))))

    # Combine treatment with fitted values
    results <- mapply(function(x, y) cbind(treatment = x, y),
      treat_by, fitted_values,
      SIMPLIFY = FALSE
    )

    # convert treatment back to original factor coding
    results <- as.data.frame(do.call(rbind, results))
    results[, "treatment"] <- factor(results[, "treatment"],
      levels = treat_levels,
      labels = treat_labels
    )
    results <- as.data.frame(results)
  } else if (isS4(fitted_object) && method == "vglm") {
    # Handle S4 object for `vglm`
    results <- as.data.frame(VGAM::fitted.values(fitted_object))
    results <- cbind(treatment = args[["treat"]], results)
  } else {
    # Default case for non-S4 object
    results <- as.data.frame(fitted_object$fitted.values)
    results <- cbind(treatment = args[["treat"]], results)
  }

  # reset rownames
  rownames(results) <- NULL

  # assign attributes and class in one go + most specific class first
  results <- structure(
    results,
    original_data  = args[[".data"]], # data actually used for fitting
    function_call  = call,
    class          = c("gps", "data.frame")
  )

  ## returning gps matrix
  return(results)
}

# internal helper: plot data.frame with a header
.print_gps_core <- function(x, ...) {
  # x: gps/csr object (data.frame with first column 'treatment')

  n <- nrow(x)
  gps_cols <- setdiff(names(x), "treatment")
  k <- length(gps_cols)

  treatment <- x[["treatment"]]
  if (is.factor(treatment)) {
    trt_levels <- levels(treatment)
  } else {
    trt_levels <- sort(unique(treatment))
  }

  cli::cli_ul()
  cli::cli_li("Number of units: {n}")
  cli::cli_li("Number of treatments: {k}")
  cli::cli_li("Treatment column: {.field treatment}")
  cli::cli_li("GPS probability columns: {.field {paste(gps_cols,
              collapse = ', ')}}")
  cli::cli_li("Treatment levels: {paste(trt_levels, collapse = ', ')}")
  cli::cli_li("All columns except 'treatment' store probabilities in [0, 1].")
  cli::cli_end()

  cli::cli_text("")

  # print underlying data.frame directly
  base::print.data.frame(x, ...)

  invisible(x)
}

#' @export
print.gps <- function(x, ...) {
  cli::cli_text("{.strong gps object} (generalized propensity scores)")
  .print_gps_core(x, ...)
}

#' @export
str.gps <- function(object, ...) {
  n <- nrow(object)
  p <- ncol(object)

  treatment <- object[["treatment"]]
  if (is.factor(treatment)) {
    trt_levels <- levels(treatment)
  } else {
    trt_levels <- sort(unique(treatment))
  }

  gps_cols <- setdiff(names(object), "treatment")
  k <- length(trt_levels)

  original_data <- attr(object, "original_data")
  call <- attr(object, "function_call")

  ## Header -------------------------------------------------------------
  cat("gps object: generalized propensity scores\n")
  cat(sprintf(" Dimensions: %d rows x %d columns\n", n, p))
  cat(sprintf(" Treatment column: %s\n", "treatment"))
  cat(sprintf(
    " GPS columns: %s\n",
    if (length(gps_cols)) paste(gps_cols, collapse = ", ") else "<none>"
  ))
  cat(sprintf(" Number of treatment levels: %d\n", k))
  cat(sprintf(
    " Treatment levels: %s\n",
    if (length(trt_levels)) paste(trt_levels, collapse = ", ") else "<none>"
  ))

  if (!is.null(original_data)) {
    cat(sprintf(
      " original_data: data.frame with %d rows and %d columns\n",
      NROW(original_data), NCOL(original_data)
    ))
  } else {
    cat(" original_data: <none>\n")
  }

  if (!is.null(call)) {
    cat(" Call:\n")
    cat("  ",
      paste(deparse(call, width.cutoff = 80L), collapse = "\n  "),
      "\n",
      sep = ""
    )
  }

  cat("\nUnderlying data.frame structure:\n")

  ## Delegate to data.frame method for the actual structure -------------
  utils::str(unclass(object), ...)

  invisible(object)
}

#' @export
summary.gps <- function(object, ...) {
  # ensure treatment is a factor
  treatment <- object[["treatment"]]
  if (!is.factor(treatment)) {
    treatment <- factor(treatment)
  }

  gps_cols <- setdiff(names(object), "treatment")

  # summary of GPS matrix by treatment level
  gps_by_trt <- lapply(
    split(object[, gps_cols, drop = FALSE], treatment),
    function(df) {
      stats_list <- lapply(df, function(v) {
        if (is.numeric(v)) {
          summary(v, ...)
        } else {
          summary(as.numeric(v), ...)
        }
      })
      stats_mat <- do.call(cbind, stats_list)
      as.matrix(stats_mat)
    }
  )

  # keep original_data as is, no summary here
  od <- attr(object, "original_data")

  res <- list(
    call             = attr(object, "function_call"),
    n                = nrow(object),
    treatments       = levels(treatment),
    gps_by_treatment = gps_by_trt,
    original_data    = od
  )
  class(res) <- "summary.gps"
  res
}

#' @export
print.summary.gps <- function(x, ...) {
  cli::cli_h1("Summary of gps object")

  # call + basic info
  if (!is.null(x$call)) {
    cli::cli_h2("Call")
    txt <- deparse(x$call, width.cutoff = 80L)
    for (line in txt) {
      cli::cli_text("{.code {line}}")
    }
    cli::cli_text("")
  }

  cli::cli_h2("Overview")
  cli::cli_ul()
  cli::cli_li("Number of units: {x$n}")
  cli::cli_li("Number of treatments: {length(x$treatments)}")
  cli::cli_li("Treatment levels: {paste(x$treatments, collapse = ', ')}")
  cli::cli_end()
  cli::cli_text("")

  ## 1) GPS matrix by treatment -----------------------------------------------
  cli::cli_h2("GPS matrix by treatment")

  gps_by_trt <- x$gps_by_treatment

  if (length(gps_by_trt) == 0L) {
    cli::cli_text("{.emph No GPS probability columns found.}")
  } else {
    for (trt in names(gps_by_trt)) {
      cli::cli_h3("Treatment: {trt}")
      stats_mat <- gps_by_trt[[trt]]

      # Row names = statistic (Min., 1st Qu., Median, Mean, 3rd Qu., Max.)
      df <- as.data.frame(stats_mat)
      df <- cbind(Statistic = rownames(stats_mat), df)
      rownames(df) <- NULL

      # round numeric columns
      num_cols <- vapply(df, is.numeric, logical(1L))
      df[num_cols] <- lapply(df[num_cols], function(v) round(v, 3))

      txt <- utils::capture.output(print(df, row.names = FALSE))
      cli::cli_verbatim(paste(txt, collapse = "\n"))
      cli::cli_text("")
    }
  }

  ## 2) original data summary --------------------------------------------------
  cli::cli_h2("Original data summary (attribute 'original_data')")

  od <- x$original_data
  # just paste the output of original summary()
  if (is.null(od)) {
    cli::cli_text("{.emph No original_data attribute or summary available.}")
  } else {
    txt <- utils::capture.output(summary(od))
    cli::cli_verbatim(paste(txt, collapse = "\n"))
  }

  invisible(x)
}

#' @export
as.data.frame.gps <- function(x, ...) {
  class(x) <- "data.frame"
  NextMethod()
}

#' @export
plot.gps <- function(x,
                     gps_cols = NULL,
                     treatment_col = "treatment",
                     ...) {
  # x: gps object (data.frame with GPS columns + treatment)
  # gps_cols: optional character vector of GPS columns to plot
  # treatment_col: name of treatment column in x
  # ...: additional arguments passed to boxplot()

  # --- basic checks ---
  trt <- x[[treatment_col]]
  if (!is.factor(trt)) {
    trt <- factor(trt)
  }

  all_gps_cols <- setdiff(names(x), treatment_col)

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

  # --- layout: number of panels based on number of treatment levels ---
  levs <- levels(trt)
  n_trt <- length(levs)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  if (n_trt == 1L) {
    graphics::par(mfrow = c(1, 1))
  } else if (n_trt == 2L) {
    graphics::par(mfrow = c(1, 2))
  } else if (n_trt <= 4L) {
    graphics::par(mfrow = c(2, 2))
  } else {
    nrow <- ceiling(n_trt / 2)
    graphics::par(mfrow = c(nrow, 2))
  }

  # --- draw boxplots for each treatment level ---
  for (lev in levs) {
    idx <- trt == lev
    dat <- x[idx, gps_cols, drop = FALSE]

    graphics::boxplot(
      dat,
      main = paste("GPS by column | treatment:", lev),
      ylab = "GPS probability",
      xlab = "GPS columns",
      ...
    )
  }

  invisible(x)
}
