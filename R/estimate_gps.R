#' Estimate Generalized Propensity Scores with various methods
#'
#' @description
#'
#' @param formula description
#' @param data description
#' @param method description
#' @param reference description
#' @param by description
#' @param missing description
#' @param fit.object description
#' @param verbose.output description
#' @param ... description
#'
#' @returns
#'
#' @details
#' Additional details...
#'
#' @examples
#' @export
estimate_gps <- function(formula,
                         data = NULL,
                         method = "multinom",
                         link = NULL,
                         reference = NULL,
                         by = NULL,
                         missing = NULL,
                         subset = NULL, # unprocessed
                         ordinal.treat = NULL, # unprocessed
                         fit.object = FALSE,
                         verbose.output = FALSE,
                         ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  call <- match.call()
  args <- list(...)

  # formula
  if (missing(formula)) {
    chk::abort_chk("The argument `formula` is missing with no default")
  }

  if (!rlang::is_formula(formula, lhs = TRUE)) {
    chk::abort_chk("The argument `formula` has to be a valid R formula with
                   treatment and predictor variables")
  }

  data.list <- .get_formula_vars(formula, data)

  args["treat"] <- list(data.list[["treat"]])
  args["covs"] <- list(data.list[["model_covs"]])

  if (is.null(args["treat"])) {
    chk::abort_chk("No treatment variable was specified")
  }

  if (is.null(args["covs"])) {
    chk::abort_chk("No predictors were specified")
  }

  if (length(args[["treat"]]) != nrow(args[["covs"]])) {
    chk::abort_chk("The treatment variable and predictors ought to have the
                   same number of samples")
  }

  if (anyNA(args[["treat"]])) {
    chk::abort_chk("The `treatment` variable can not have any NA's")
  }

  n_levels <- nunique(args[["treat"]])
  if (n_levels > 10) {
    chk::wrn("The `treatment` variable has more than 10 unique levels. Consider
             dropping the number of groups, as the vector matching algorithm may
             not perform well")
  }

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
    args[['treat']] <- factor(args[['treat']], levels = unique(args[['treat']], ordered = FALSE))
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

  if(!is.null(args[['link']]) || !missing(args[['link']])) {
    chk::chk_string(args[['link']])

    if(args[['link']] %nin% available_links) {
      chk::abort_chk(sprintf('The argument `link` for the method %s only accepts values: %s',
                             method, word_list(add_quotes(available_links))))
    }
  } else {
    args[['link']] <- available_links[1]
  }

  # reference
  levels_treat <- as.character(unique(args[["treat"]]))

  if(!is.null(ordinal.treat) && !is.null(reference)) {
    chk::wrn('There is no need to specify `reference` if `ordinal.treat` was provided. Ignoring the `reference` argument')
  } else {
    if (is.null(reference)) {
      reference <- levels_treat[1]
      args[['treat']] <- stats::relevel(args[['treat']], ref = reference)
    } else if (!(is.character(reference) && length(reference) == 1L && !anyNA(reference))) {
      chk::abort_chk("The argument `reference` must be a single string of length 1")
    } else if (!(reference %in% levels_treat)) {
      chk::abort_chk("The argument `reference` is not in the unique levels of the
                   treatment variable")
    } else {
      args[['treat']] <- stats::relevel(args[['treat']], ref = reference)
    }
  }

  # missing
  missing <- .process_missing(missing, method)

  # by --> assembling the arguments

  # fit.object + verbose.output
  chk::chk_all(list(fit.object, verbose.output), chk::chk_flag)

  # assembling the arguments list
  args["covs"] <- list(data.list[["reported_covs"]])
  args[".formula"] <- list(formula)
  args[".data"] <- list(data)
  args["method"] <- list(method)
  args["reference"] <- reference
  args[["missing"]] <-
  args[["by"]] <- .process_by(by, data, args[["treat"]])
  args["fit.object"] <- list(fit.object)
  args["verbose.output"] <- list(verbose.output)
  args["subset"] <- list(subset)

  ####################### FITTING ##############################################
  fit.func <- .gps_methods[[method]]$func_used

  fitted_object <- do.call(
    fit.func,
    args
  )

  ####################### OUTPUT OBJECT ########################################
  # defining the class of the output
  if (fit.object) {
    return(fitted_object)
  } else if (isS4(fitted_object)) {
    if (method == "vglm") {
      gps <- VGAM::fitted.values(fitted_object)
      class(gps) <- "gps"
      return(gps)
    }
  } else {
    gps <- as.matrix(fitted_object$fitted.values)
    class(gps) <- "gps"
    return(gps)
  }
}
