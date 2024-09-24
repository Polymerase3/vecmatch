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
#'
#' @export
estimate_gps <- function(formula,
                         data = NULL,
                         method = 'multinom',
                         reference = NULL,
                         by = NULL,             ##not processed
                         missing = NULL,        ##not processed
                         fit.object = FALSE,
                         verbose.output = FALSE,
                         ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  call <- match.call()
  args <- list(...)

  #formula
  if(missing(formula)) {
    chk::abort_chk('The argument `formula` is missing with no default')
  }

  if(!rlang::is_formula(formula, lhs = TRUE)) {
    chk::abort_chk('The argument `formula` has to be a valid R formula with
                   treatment and predictor variables')
  }

  data.list <- .get_formula_vars(formula, data)
  args['treat'] <- list(data.list[['treat']])
  args['covs'] <- list(data.list[['model_covs']])
  args['.covs'] <- list(data.list[['reported_covs']])
  args['.formula'] <- list(formula)

  if(is.null(args['treat'])) {
    chk::abort_chk('No treatment variable was specified')
  }

  if(is.null(args['covs'])) {
    chk::abort_chk('No predictors were specified')
  }

  if(length(args[['treat']]) != nrow(args[['covs']])) {
    chk::abort_chk('The treatment variable and predictors ought to have the
                   same number of samples')
  }

  if(anyNA(args[['treat']])) {
    chk::abort_chk("The `treatment` variable can not have any NA's")
  }

  n_levels <- nunique(args[['treat']])
  if (n_levels > 7) {
    chk::wrn('The `treatment` variable has more than 7 unique levels. Consider
             dropping the number of groups, as the vector matching algorithm may
             not perform well')
  }

  #data
  if(!is.null(data)) .check_df(data)
  args['.data'] <- list(data)

  #method
  if(is.null(substitute(method)) || missing(method)) {
    method <- 'multinom'
    attr(method, "name") <- method
  } else {
    .check_method(method)
    method.name <- deparse1(substitute(method))
    attr(method, "name") <- method.name
  }

  args['method'] <- list(method)

  #reference
  levels_treat <- as.character(unique(args[['treat']]))

  if(is.null(reference)) {
    reference <- levels_treat[1]
  } else if (!(is.character(reference) && length(reference) == 1L && !anyNA(reference))) {
    chk::abort_chk('The argument `reference` must be a single string of length 1')
  }

  if(!(reference %in% levels_treat)) {
    chk::abort_chk('The argument `reference` is not in the unique levels of the
                   treatment variable')
  }

  args['reference'] <- reference

  # fit.object + verbose.output
  chk::chk_all(list(fit.object, verbose.output), chk::chk_logical)
  args['fit.object'] <- list(fit.object)
  args['verbose.output'] <- list(verbose.output)

  ####################### FITTING ##############################################

  ####################### OUTPUT OBJECT ########################################

  #defining the class of the output
}
