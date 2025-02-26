#' @title Calculate treatment allocation probabilities
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
#'   levels specified in the argument. Is a call to the function `factor(treat,
#'   levels = ordinal_treat, ordered = TRUE`.
#' @param fit_object a logical flag. If `TRUE`, the the fitted object is
#'   returned instead of the GPS matrix.
#' @param verbose_output a logical flag. If `TRUE` a more verbose version of the
#'   function is run and the output is printed out to the console.
#' @param ... additional arguments, that can be passed to the fitting function
#'   and are not controlled by the above arguments. For more details and
#'   examples refer to the Details section and documentations of corresponding
#'   functions.
#'
#' @return A numeric matrix of class `gps` with the number of columns equal to
#'   the number of unique treatment variable levels plus one (for the treatment
#'   variable itself) and the number of row equal to the number of subjects in
#'   the initial dataset. The original dataset used for estimation can be
#'   accessed as `original_data` attribute.
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
#' library("brglm2")
#'
#' # Conducting covariate balancing on the `airquality` dataset. Our goal was to
#' # compare ozone levels by month, but we discovered that ozone levels are
#' # strongly correlated with wind intensity (measured in mph), and the average
#' # wind intensity varies across months. Therefore, we need to balance the
#' # months by wind values to ensure a valid comparison of ozone levels.
#'
#' # Initial imbalance of means
#' tapply(airquality$Wind, airquality$Month, mean)
#'
#' # Formula definition
#' formula_air <- formula(Month ~ Wind)
#'
#' # Estimating the generalized propensity scores using brglm2 method using
#' # maximum penalized likelihood estimators with powers of the Jeffreys
#' gp_scores <- estimate_gps(formula_air,
#'   data = airquality, method = "brglm2",
#'   reference = "5", verbose_output = TRUE,
#'   control = brglmControl(type = "MPL_Jeffreys")
#' )
#'
#' # Filtering the observations outside the csr region
#' gps_csr <- csregion(gp_scores)
#'
#' # Calculating imbalance after csr
#' filter_which <- attr(gps_csr, "filter_vector")
#' filtered_air <- airquality[filter_which, ]
#'
#' tapply(filtered_air$Wind, filtered_air$Month, mean)
#'
#' # We can also investigate the imbalance using the raincloud function
#' raincloud(filtered_air,
#'   y = Wind,
#'   group = Month,
#'   significance = "t_test"
#' )
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
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  call <- match.call()
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

      # model the data
      fit <- do.call(
        fit_func,
        as.list(by_env)
      )

      # append to existing list
      fitted_object <- append(fitted_object, list(fit))

      # delete env
      rm(by_env)

      ## save treatment var
      treat_by[[i]] <- args[["treat"]][by_sub]
    }
  } else {
    fitted_object <- do.call(
      fit_func,
      args
    )
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

  # Assign attributes and class
  class(results) <- c("data.frame", "gps")

  ## adding original data as attribute
  attr(results, "original_data") <- args[[".data"]]

  ## returning gps matrix
  return(results)
}
