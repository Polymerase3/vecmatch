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
#'   Doesn't need to be provided, as the columns can also be specified directly
#'   and explicitly in the formula, using `$` operator and the name of the
#'   dataset.
#' @param method a single string describing the model used for the calculation
#'   of generalized propensity scores. The default value is set to `multinom`.
#'   For available methods refer to the Details section below
#' @param link a single string; determines an alternative model for a method
#'   used for estimation. For available links, see XXX --> Actually do wywalenia
#'   i zastąpić to method
#' @param subset a logical atomic vector of length equal to the number of rows
#'   in the `data` arguments. Allows to filter out observations from the further
#'   analysis, for which the value of the vector is qual to `FALSE`.
#' @param reference a single string describing one class from the treatment
#'   variable, referred to as the baseline category in the calculation of
#'   generalized propensity scores.
#' @param by a single string with the name of a column, contained in the `data`
#'   argument. Tha dataset will be divided by the groups created by the grouping
#'   `by` variable and the calculation of the propensity scores will be carried
#'   out separately for each group. The results will then be merged and
#'   presented to the uesr as a single GPS matrix.
#' @param missing a single string describing the mtehod to be used in the
#'   handling of the missing data. For possible values refer to XXX. The default
#'   value is `complete.cases`. FIXXXXXX
#' @param ordinal.treat an atomic vector of the length equal to the length of
#'   unique levels of the treatment variable. Confirms, that the treatment
#'   variable is an ordinal variable and adjusts its levels, to the order of
#'   levels specified in the argument. Is a call to the function `factor(treat,
#'   levels = ordinal.treat, ordered = TRUE`.
#' @param fit.object a logical flag. If `TRUE`, the the fitted object is
#'   returned instead of the GPS matrix.
#' @param verbose.output a logical flag. If `TRUE` a more verbose version of the
#'   function is run and the output is printed out to the console.
#' @param ... additional arguments, that can be passed to the fitting function
#'   and are not controlled by the above arguments. For more details and
#'   examples refer to XXX
#'
#' @returns A numeric matrix with the number of columns equal to the number of
#'   unique treatment variable levels plus one (for the treatment variable
#'   itself) and the number of row equal to the number of subjects in the
#'   initial dataset.
#'
#' @details The main goal of the `estimate_gps()` function is to calculate the
#'   generalized propensity scores aka. treatment allocation probabilities. It
#'   is the first step in the workflow vector matching algorithm and is
#'   essential for the further analysis. The returned matrix of class `gps` can
#'   then be passed to the `csregion()` function to calculate the rectangular
#'   common support region boundaries and drop samples uneligible for the
#'   further analysis. The list of available methods operated by the
#'   `estimate_gps()` is provided below with a short description and function
#'   used for the calculations:
#'   * `multinom` - multinomial logistic regression model [nnet::multinom()]
#'   * `vglm` - vector generalized linear model for multinomial data [VGAM::vglm()],
#'   * `brglm2` - bias reduction model for multinomial respones using the poisson trick [brglm2::brmultinom()],
#'   * `mblogit` - baseline-category logit models [mclogit::mblogit()].
#'   * `polr` - ordered logistic or probit regression onyl for ordered factor variables from [MASS::polr()]. The `method` argument of the underlying `MASS::polr()` package function can be controlled with the `link` argument. Available options: `link = c("logistic", "probit", "loglog", "cloglog", "cauchit")`
#' @examples
#'
#'library('brglm2')
#'
#'# Conducting covariate balancing on the `airquality` dataset. Our goal was to
#'# compare ozone levels by month, but we discovered that ozone levels are strongly
#'# correlated with wind intensity (measured in mph), and the average wind intensity
#'# varies across months. Therefore, we need to balance the months by wind values
#'# to ensure a valid comparison of ozone levels.
#'
#'# Initial imbalance of means
#'tapply(airquality$Wind, airquality$Month, mean)
#'
#'# Formula definition
#'formula_air <- formula(Month ~ Wind)
#'
#'# Estimating the generalized propensity scores using brglm2 method using
#'# maximum penalized likelihood estimators with powers of the Jeffreys
#'gp_scores <- estimate_gps(formula_air, data = airquality, method = 'brglm2',
#'                          reference = '5', verbose.output = TRUE,
#'                          control = brglmControl(type = 'MPL_Jeffreys'))
#'
#' # Filtering the observations outside the csr region
#' gps_csr <- csregion(gp_scores)
#'
#' # Calculating imbalance after csr
#' filter_which <- attr(gps_csr, 'filter_vector')
#' filtered_air <- airquality[filter_which, ]
#'
#' tapply(filtered_air$Wind, filtered_air$Month, mean)
#'
#' # We can also investigate the imbalance using the raincloud function
#' raincloud(filtered_air,
#'           y = Wind,
#'           group = Month,
#'           significance = 't_test')
#' @export

estimate_gps <- function(formula,
                         data = NULL,
                         method = "multinom",
                         link = NULL,
                         reference = NULL,
                         by = NULL,
                         missing = NULL, # unprocessed
                         subset = NULL,
                         ordinal.treat = NULL, # unprocessed
                         fit.object = FALSE,
                         verbose.output = FALSE,
                         ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  call <- match.call()
  args <- list(...)

  dots <- substitute(list(...))[-1]

  ## If function call then substitute, else normal
  additional_args <- which(names(call)[-1] %nin% names(formals(estimate_gps)))

  if(!is.null(dots)) {
    for(i in seq_along(additional_args)) {
      argname <- names(args)[i]
      callname <- paste0(argname, '_call')
      callname_char <- paste0(callname, '_char')
      args[callname] <- FALSE

      if(is.call(dots[[i]])) {
        args[callname] <- TRUE
        args[callname_char] <- as.character(dots[[i]])[1]
      }
    }
  }

  # formula
  data.list <- .process_formula(formula, data)

  # args assignment to list used in calculations
  args["treat"] <- list(data.list[["treat"]])
  args["covs"] <- list(data.list[["model_covs"]])


  # process and check ordinal.treat
  if (!is.null(ordinal.treat)) {
    chk::chk_atomic(ordinal.treat)
    chk::chk_vector(ordinal.treat)

    if (length(ordinal.treat) != length(unique(args[["treat"]]))) {
      chk::abort_chk("The numbers of levels provided in `ordinal.treat` has to
                     be the same, as the number of unique levels in the treatment
                     variable.")
    }
    args[["treat"]] <- factor(args[["treat"]], levels = ordinal.treat, ordered = TRUE)
  } else {
    args[["treat"]] <- factor(args[["treat"]], levels = unique(args[["treat"]], ordered = FALSE))
  }

  # data
  if (!is.null(data)) .check_df(data)

  # method
  if (is.null(substitute(method)) || missing(method)) {
    method <- "multinom"
    attr(method, "name") <- method
  } else {
    .check_method(method)
    method.name <- deparse1(substitute(method))
    attr(method, "name") <- method.name
  }

  # link
  available_links <- .gps_methods[[method]]$link_fun

  if (!is.null(link)) {
    chk::chk_string(link)

    if (link %nin% available_links) {
      chk::abort_chk(sprintf(
        "The argument `link` for the method %s only accepts values: %s",
        method, word_list(add_quotes(available_links))
      ))
    }

    args[["link"]] <- link
  } else {
    args[["link"]] <- available_links[1]
  }

  # reference
  ref.list <- .process_ref(args[['treat']],
                           ordinal.treat = ordinal.treat,
                           reference = reference)

  args[['treat']] <- ref.list[['data.relevel']]
  reference <- ref.list[['reference']]

  # missing
  missing <- .process_missing(missing, method)

  # subset
  if (!is.null(subset)) {
    chk::chk_string(subset)

    if (subset %nin% colnames(data)) {
      chk::abort_chk(sprintf(
        "The column %s defined in the `subset` argument was not found in the provided dataset.",
        subset
      ))
    }

    subset_logvec <- as.vector(data[[subset]])
    if (!is.logical(subset_logvec) || length(dim(subset_logvec)) == 2L) {
      chk::abort_chk("The `subset` argument has to be a name of single column with logical values.")
    }

    use.subset <- TRUE
  } else {
    use.subset <- FALSE
  }

  # fit.object + verbose.output
  chk::chk_all(list(fit.object, verbose.output), chk::chk_flag)

  # assembling the arguments list
  if (use.subset) {
    args[["treat"]] <- args[["treat"]][subset_logvec]
    args["covs"] <- list(data.list[["reported_covs"]][subset_logvec, ])

    args[".data"] <- list(data[subset_logvec, ])
    args[["by"]] <- .process_by(by, data, args[["treat"]])[subset_logvec, ]
  } else {
    args["covs"] <- list(data.list[["reported_covs"]])
    args[".data"] <- list(data)
    args[["by"]] <- .process_by(by, data, args[["treat"]])
  }

  args["formula"] <- list(formula)
  args["method"] <- list(method)
  args["reference"] <- reference
  args[["missing"]] <- missing
  args["fit.object"] <- list(fit.object)
  args["verbose.output"] <- list(verbose.output)
  # args["subset"] <- list(subset)

  ####################### FITTING ##############################################
  fit.func <- .gps_methods[[method]]$func_used
  use.by <- FALSE
  if (!is.null(args[["by"]])) {
    fitted_object <- list()
    by.levels <- levels(attr(args[["by"]], "by.factor"))
    use.by <- TRUE

    for (i in seq_along(by.levels)) {
      # subset rule
      by.sub <- attr(args[["by"]], "by.factor") == by.levels[i]

      # create env and subset vars
      by.env <- list2env(args, envir = new.env(), parent = emptyenv())
      with(by.env, {
        selected <- mget(c(".data", "covs", "treat"), envir = by.env)
        subsetted <- lapply(selected, function(x) {
          if (is.data.frame(x)) {
            x[by.sub, ]
          } else if (is.atomic(x)) {
            x[by.sub]
          }
        })
      })

      # overwrite
      list2env(by.env$subsetted, envir = by.env)

      # model the data
      fit <- do.call(
        fit.func,
        as.list(by.env)
      )

      # append to existing list
      fitted_object <- append(fitted_object, list(fit))

      # delete env
      rm(by.env)
    }

  } else {
    fitted_object <- do.call(
      fit.func,
      args
    )
  }

  ####################### OUTPUT OBJECT ########################################
  # defining the class of the output
  if (fit.object) {
    return(fitted_object)
  } else if (use.by) {
    if(all(vapply(fitted_object, isS4, logical(1L)))) {
      if(method == 'vglm') {
        results <- lapply(fitted_object, VGAM::fitted.values)
        gps <- do.call(rbind, results)
      }
    } else {
      results <- lapply(fitted_object, '[[', 'fitted.values')
      gps <- do.call(rbind, results)
    }
  } else if (isS4(fitted_object)) {
    if (method == "vglm") {
      gps <- VGAM::fitted.values(fitted_object)
    }
  } else {
    gps <- fitted_object$fitted.values
    gps <- as.data.frame(cbind(treatment = args[['treat']], gps))
  }

  ## cbind treatment to all outputs!
  ## add argnames as attributes!

  class(gps) <- c('data.frame', 'gps')
  return(gps)
}
