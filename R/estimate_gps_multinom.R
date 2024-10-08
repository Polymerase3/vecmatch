.estimate_gps_multinom <- function(.formula, treat, covs, method,
                                   fit.object, verbose.output,
                                   subset, ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  fit_object <- NULL
  Args <- list(...)


  if (!is.null(subset)) {
    covs <- covs[subset, , drop = FALSE]
    treat <- treat[subset]
  }

  ## Assign and check treatment type
  treat <- .assign_treatment_type(treat)
  treat.type <- .get_treat_type(treat)

  ## Process ordinal treat and scale

  data <- data.frame(treat, covs)
  data <- as.data.frame(lapply(data, scale_0_to_1))

  ## Defining the list of arguments and processing
  Args[['formula']] <- .formula
  Args[['treat']] <- treat
  Args[['covs']] <- covs
  Args[['data']] <- data
  Args[['method']] <- method
  Args[['model']] <- fit.object
  Args[['verbose']] <- verbose.output

  ## Fit the model
  if(method == 'multinom') {
    infos <- .gps_methods[['multinom']]
    rlang::check_installed(infos$packages_needed)

    ## Processing the stuff
    Args <- match_add_args(arglist = Args,
                           funlist = infos$fun.arg.check) ## check if works cause funs.arg.check quoted







  }
}
