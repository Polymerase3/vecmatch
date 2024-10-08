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
                         method = 'glm',
                         reference = NULL,
                         by = NULL,
                         missing = NULL,
                         subset = NULL,     #unprocessed
                         ordinal.treat = NULL,    #unprocessed
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

  #process and check ordinal.treat
  if(!is.null(ordinal.treat)) {
    chk::chk_atomic(ordinal.treat)
    chk::chk_vector(ordinal.treat)

    if(length(ordinal.treat) != length(unique(args[['treat']]))) {
      chk::abort_chk('The numbers of levels provided in `ordinal.treat` has to
                     be the same, as the number of unique levels in the treatment
                     variable.')
    }
    args[['treat']] <- factor(args[['treat']], levels = ordinal.treat, ordered = TRUE)
  }

  if(is.null(ordinal.treat)) {
    args[['treat']] <- factor(args[['treat']], ordered = FALSE)
  }

  #data
  if(!is.null(data)) .check_df(data)

  #method
  if(is.null(substitute(method)) || missing(method)) {
    method <- 'multinom'
    attr(method, "name") <- method
  } else {
    .check_method(method)
    method.name <- deparse1(substitute(method))
    attr(method, "name") <- method.name
  }

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

  #missing
  if(is.null(missing)) {
    missing <- 'complete.cases'
  } else if (!(is.character(missing) && length(missing) == 1L && !anyNA(missing))) {
    chk::abort_chk('The argument `missing` must be a single string of length 1')
  }

  #by --> assembling the arguments

  # fit.object + verbose.output
  chk::chk_all(list(fit.object, verbose.output), chk::chk_logical)

  # assembling the arguments list
  args['covs'] <- list(data.list[['reported_covs']])
  args['.formula'] <- list(formula)
  args['.data'] <- list(data)
  args['method'] <- list(method)
  args['reference'] <- reference
  args[['missing']] <- .process_missing(missing, method)
  args[['by']] <- .process_by(by, data, args[['treat']])
  args['fit.object'] <- list(fit.object)
  args['verbose.output'] <- list(verbose.output)
  args['subset'] <- list(subset)

  ####################### FITTING ##############################################
  fit.func <- .gps_methods[[method]]$func_used

  fitted_object <- do.call(fit.func,
                           args)

  ####################### OUTPUT OBJECT ########################################
  #defining the class of the output
  if(fit.object) {
    return(fitted_object)
  } else {
    gps <- as.matrix(fitted_object$fitted.values)
    class(gps) <- 'gps'
    return(gps)
  }
}
