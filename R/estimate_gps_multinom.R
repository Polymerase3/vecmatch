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
  Args[["model"]] <- fit.object
  Args[["verbose"]] <- verbose.output

  ## Fit the model
  if (treat.type == "multinom" || treat.type == "binary") {
    if (method == "multinom") {
      infos <- .gps_methods[["multinom"]]
      rlang::check_installed(infos$packages_needed)
      sapply(infos$packages_needed, requireNamespace, quietly = TRUE)

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

    if (method == "vglm") {
      infos <- .gps_methods[["vglm"]]
      rlang::check_installed(infos$packages_needed)
      sapply(infos$packages_needed, requireNamespace, quietly = TRUE)

      ## Check link
      if (is.null(Args[["link_fun"]])) {

        link_fun <- infos$link_fun[1]
        chk::message_chk(sprintf("You can specify the type of %s model using the `link_fun` argument. \n
                                 The default value %s was set.", method, Args[["link_fun"]]))
      } else {
        if (Args[["link_fun"]] %nin% infos$link_fun) {
          chk::abort_chk(sprintf(
            "The argument `link_fun` only accepts following values: %s",
            paste(infos$link_fun, collapse = ", ")
          ))
        } else {
          link_fun <- Args[['link_fun']]
        }
      }

      ## Processing Args
      Args <- match_add_args(
        arglist = Args,
        infos$fun.arg.check
      )

      ## Overwriting Args
      Args[['family']] <- VGAM::multinomial
      Args[['trace']] <- verbose.output
      Args[['control']] <- VGAM::vglm.control(...)

      ## Fit model
      if (link_fun %in% infos$link_fun) {
        tryCatch(
          verbosely(
            {
              fit <- do.call(switch(link_fun,
                                    'multinomial_logit' = VGAM::vglm,
                                    'reduced_rank_ml' = VGAM::rrvglm),
                args = Args
              )
            },
            verbose = verbose.output
          ),
          error = function(e) {
            chk::abort_chk(sprintf(
              "There was a problem fitting the %s regressions with `VGAM::vglm(family = multinom)`.\n
                               Error message: (from `VGAM::vglm(family = multinom)`) %s",
              link_fun, conditionMessage(e)
            ), tidy = FALSE)
          }
        )

        fitted_obj <- fit
      } else {
        probably_a_bug <- TRUE
      }
    }
  } else {
    fitted_obj <- NULL
  }

  if (!is.null(fitted_obj)) {
    return(fitted_obj)
  } else {
    probably_a_bug <- TRUE
  }

  if (probably_a_bug) {
    chk::abort_chk("The function `estimate_gps()` was not able to estimate the propensity scores.
                   It's probably a bug. Please let the author know.")
  }
}
