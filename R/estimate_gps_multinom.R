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

  ## Process method (if method and treat type)
  # if ()


  ## Process data
  data <- data.frame(treat, covs)
  data <- as.data.frame(lapply(data, scale_0_to_1))

  ## Defining the list of arguments and processing
  Args[["formula"]] <- update.formula(.formula, treat ~ .)
  Args[["treat"]] <- treat
  Args[["covs"]] <- covs
  Args[["data"]] <- data
  Args[["method"]] <- method
  Args[["model"]] <- fit.object
  Args[["verbose"]] <- verbose.output

  ## Fit the model
  if (treat.type == "multinom" || treat.type == "binary") {
    if (method == "multinom") {
      infos <- .gps_methods[["multinom"]]
      rlang::check_installed(infos$packages_needed)

      ## Processing the stuff
      Args <- match_add_args(
        arglist = Args,
        funlist = infos$fun.arg.check
      )

      ## Fit the multinom
      tryCatch(
        verbosely(
          {
            fit <- do.call(nnet::multinom,
              args = Args,
              quote = TRUE
            )
          },
          verbose = verbose.output
        ),
        error = function(e) {
          chk::abort_chk(sprintf(
            "There was a problem fitting the multinomial %s regressions with `nnet::multinom()`.\n
                               Error message: (from `nnet::multinom()`) %s",
            infos[["link_fun"]], conditionMessage(e)
          ), tidy = FALSE)
        }
      )

      fitted_obj <- fit
    }
  } else {
    fitted_obj <- NULL
  }

  if (!is.null(fitted_obj)) {
    return(fitted_obj)
  } else {
    chk::abort_chk("The function `estimate_gps()` was not able to estimate the propensity scores.
                   It's probably a bug. Please let the author know.")
  }
}
