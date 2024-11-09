.estimate_gps_multinom <- function(formula, treat, link, covs, method,
                                   fit.object, verbose.output,
                                   subset, ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  fit_object <- NULL
  Args <- list(...)
  probably_a_bug <- FALSE

  ## Assign and check treatment type
  treat <- .assign_treatment_type(treat)
  treat.type <- .get_treat_type(treat)

  ## Check polr with treatment type
  if(method == 'polr' && treat.type != 'ordinal') {
    chk::abort_chk('If `method = "polr"`, the treatment variable must be an ordered factor.')
  }

  ## Process data
  data <- data.frame(treat, covs)
  data <- as.data.frame(lapply(data, scale_0_to_1))

  ## Defining the list of arguments and processing
  if (is.atomic(covs)) {
    Args[["formula"]] <- stats::update.formula(formula, treat ~ covs)
  } else {
    Args[["formula"]] <- stats::update.formula(formula, treat ~ .)
  }
  Args[["treat"]] <- treat
  Args[["covs"]] <- covs
  Args[["data"]] <- data
  Args[["model"]] <- fit.object
  Args[["verbose"]] <- verbose.output
  Args[['link']] <- link

  ################### FITTING THE MODELS #######################################
  if (treat.type == "multinom" || treat.type == "binary" || (treat.type == "ordinal" && method != 'polr')) {
    ## --NNET::multinom()--------------------------------------------------------
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
            infos[["link"]], conditionMessage(e)
          ), tidy = FALSE)
        }
      )

      fitted_obj <- fit
    }

    ## --VGLM--------------------------------------------------------------------
    if (method == "vglm") {
      infos <- .gps_methods[["vglm"]]
      rlang::check_installed(infos$packages_needed)
      sapply(infos$packages_needed, requireNamespace, quietly = TRUE)

      ## Procesing additional args
      ## Control
       if(is.null(Args[['control']])) {
         if(link == "multinomial_logit") {
           Args[["control"]] <- VGAM::vglm.control()
         } else if (link == 'reduced_rank_ml') {
           Args[["control"]] <- VGAM::rrvglm.control()
         }

       } else {
         if(link == "multinomial_logit") {
           if(!Args[['control_call']] ||
              (Args[['control_call_char']] != "VGAM::vglm.control" &&
               Args[['control_call_char']] != 'vglm.control')) {
             chk::abort_chk('The argument control has to be a valid function call to the function
                          VGAM::vglm.control()')
           }
         } else if(link == 'reduced_rank_ml') {
           if(!Args[['control_call']] ||
              (Args[['control_call_char']] != "VGAM::rrvglm.control" &&
               Args[['control_call_char']] != 'rrvglm.control')) {
             chk::abort_chk('The argument control has to be a valid function call to the function
                          VGAM::rrvglm.control()')
           }
         }
       }

      ## family
      if(is.null(Args[['family']])) {
        Args[['family']] <- VGAM::multinomial()
      } else {
        tryCatch({
          VGAM::vglm(Args[['formula']], family = Args[['family']], data = Args[['data']])
        }, error = function(e) {
          chk::abort_chk('The `family` argument has to be a valid VGAM family function argument.')
        })
      }

      ## Processing Args
      Args <- match_add_args(
        arglist = Args,
        infos$fun.arg.check
      )

      ## Overwriting Args
      ## trace (verose)
      Args[["trace"]] <- verbose.output

      ## method
      if(link == 'reduced_rank_ml') Args[['method']] <- 'rrvglm.fit'

      fun_used <- ifelse(link == "multinomial_logit", "`VGAM::vglm()`",
        "`VGAM::rrvglm()`"
      )

      ## Fit model
      if (link %in% infos$link_fun) {
        tryCatch(
          verbosely(
            {
              fit <- do.call(
                switch(link,
                  "multinomial_logit" = VGAM::vglm,
                  "reduced_rank_ml" = VGAM::rrvglm
                ),
                args = Args
              )
            },
            verbose = verbose.output
          ),
          error = function(e) {
            chk::abort_chk(sprintf(
              "There was a problem fitting the %s regressions with %s.\n
                               Error message: (from %s) %s",
              link, fun_used, fun_used, conditionMessage(e)
            ), tidy = FALSE)
          }
        )

        fitted_obj <- fit
      } else {
        probably_a_bug <- TRUE
      }
    }

    ## --brglm2::brmultinom()------------------------------------------------------------
    ####################################### TO WORK ON########################################
    if (method == "brglm2") {
      infos <- .gps_methods[["brglm2"]]
      rlang::check_installed(infos$packages_needed)
      sapply(infos$packages_needed, requireNamespace, quietly = TRUE)

      ## Processin control arg
      if(is.null(Args[['control']])) {
        Args[['control']] <- brglm2::brglmControl()
      } else {
        if(link == "baseline_category_logit") {
          if(!Args[['control_call']] ||
             (Args[['control_call_char']] != "brglm2::brglmControl" &&
              Args[['control_call_char']] != 'brglmControl')) {
            chk::abort_chk('The argument control has to be a valid function call to the function
                          brglm2::brglmControl()')
          }
        }
      }

      ## Processing the arguments
      Args <- match_add_args(
        arglist = Args,
        infos$fun.arg.check
      )

      ## Fit the brglm2
      if(link %in% infos$link_fun) {
        tryCatch(
          verbosely(
            {
              fit <- do.call(brglm2::brmultinom,
                             args = Args
              )
            },
            verbose = verbose.output
          ),
          error = function(e) {
            chk::abort_chk(sprintf(
              "There was a problem fitting the %s regressions with `brglm2::brmultinom()`.\n
                               Error message: (from `brglm2::brmultinom()`) %s",
              link, conditionMessage(e)
            ), tidy = FALSE)
          }
        )

        fitted_obj <- fit
      } else {
        probably_a_bug <- TRUE
      }
    }

    ## --mclogit::mclogit()------------------------------------------------------------
    if (method == "mblogit") {
      infos <- .gps_methods[["mblogit"]]
      rlang::check_installed(infos$packages_needed)
      sapply(infos$packages_needed, requireNamespace, quitely = TRUE)
      names.matrix <- unique(Args[['treat']])
      ncol.matrix <- length(names.matrix)

      ## Process the control arg
      if(is.null(Args[['control']])) {
        Args[['control']] <- mclogit::mclogit.control()
      } else {
        if(link == "baseline_category_logit" || "conditional_logit") {
          if(!Args[['control_call']] ||
             (Args[['control_call_char']] != "mclogit::mclogit.control" &&
              Args[['control_call_char']] != 'mclogit.control' &&
              Args[['control_call_char']] != "mclogit::mmclogit.control" &&
              Args[['control_call_char']] != 'mmclogit.control')
             ) {
            chk::abort_chk('The argument control has to be a valid function call to the function
                          mclogit::mclogit.control()')
          }
        }
      }

      ## Overwriting Args
      if (is.null(Args[["estimator"]])) {
        Args[["estimator"]] <- "ML"
      }

      ## Processing Args
      Args <- match_add_args(
        arglist = Args,
        infos$fun.arg.check
      )

      if(link %in% infos$link_fun) {
        # Fit model
        tryCatch(
          verbosely(
            {
              fit <- do.call(mclogit::mblogit,
                             args = Args
              )
            },
            verbose = verbose.output
          ),
          error = function(e) {
            chk::abort_chk(sprintf(
              "There was a problem fitting the %s regressions with %s.\n
                               Error message: (from %s) %s",
              link, fun_used, fun_used, conditionMessage(e)
            ), tidy = FALSE)
          }
        )

        fitted_obj <- fit

        ## processing the fitted.values matrix
        fitted.matrix <- matrix(NA, nrow = length(fitted_obj$fitted.values)/ncol.matrix,
                                ncol = ncol.matrix)
        colnames(fitted.matrix) <- names.matrix
        for (i in 1:ncol.matrix) {
          sub_vector <- seq(from = i, to = length(fitted_obj$fitted.values), by = ncol.matrix)

          fitted.matrix[, i] <- fitted_obj$fitted.values[sub_vector]
        }

        fitted_obj$fitted.values <- fitted.matrix

      } else {
        probably_a_bug <- TRUE
      }
    }

    ## --thats where ordinals should start
  } else if (treat.type == "ordinal") {
    if(method == 'polr') {
      infos <- .gps_methods[['polr']]
      rlang::check_installed(infos$packages_needed)

      # map link to method arg of MASS::polr
      which_change <- which(names(Args) == 'link')
      names(Args)[which_change] <- 'method'


      ## Processing Args
      Args <- match_add_args(
        arglist = Args,
        infos$fun.arg.check
      )

      if(link %in% infos$link_fun) {
        tryCatch(
          verbosely(
            {
              fit <- do.call(MASS::polr,
                             args = Args
              )
            },
            verbose = verbose.output
          ),
          error = function(e) {
            chk::abort_chk(sprintf(
              "There was a problem fitting the %s regressions with `MASS::polr()`.\n
                               Error message: (from `MASS::polr()`) %s",
              link, conditionMessage(e)
            ), tidy = FALSE)
          }
        )
      }
        fitted_obj <- fit
    }
  } else {
    probably_a_bug <- TRUE
  }

  ## Last check of output object
  if (probably_a_bug) {
    chk::abort_chk("The function `estimate_gps()` was not able to estimate the propensity scores.
                   It's probably a bug. Please let the author know.")
  } else {
    return(fitted_obj)
  }
}
