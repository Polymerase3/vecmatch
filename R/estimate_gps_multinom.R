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

  data <- data.frame(factor(treat), covs)
  data <- as.data.frame(lapply(data, scale_0_to_1))

  ## Defining the list of arguments and processing
  Args[['formula']] <- .formula
  Args[['treat']] <- treat
  Args[['covs']] <- covs
  Args[['data']] <- data
  Args[['method']] <- method
  Args[['model']] <- fit.object
  Args[['verbose']] <- verbose.output

  ## Processing the stuff
  Args <- match_add_args(arglist = Args,
                         funlist = list(nnet::multinom,
                                        nnet::nnet.default,
                                        nnet::nnet.formula))

  ## Fit the model
  if(method == 'multinom') {
    .gps_methods[['multinom']]
    rlang::check_installed("nnet")



  }
}
