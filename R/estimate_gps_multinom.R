.estimate_gps_multinom <- function(formula, treat, link, covs, method,
                                   fit_object, verbose_output,
                                   subset, ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  fit_object <- NULL
  args <- list(...)
  probably_a_bug <- FALSE

  ## Assign and check treatment type
  treat <- .assign_treatment_type(treat)
  treat_type <- .get_treat_type(treat)

  ## Check polr with treatment type
  .chk_cond(
    method == "polr" && treat_type != "ordinal",
    'If `method = "polr"`, the treatment variable
            must be an ordered factor. Use the `ordinal_treat` argument to
    define levels.'
  )

  ## Process data
  covs <- as.data.frame(lapply(covs, scale_0_to_1))
  data <- data.frame(treat, covs)

  ## Defining the list of arguments and processing
  if (is.atomic(covs)) {
    args[["formula"]] <- stats::update.formula(formula, treat ~ covs)
  } else {
    args[["formula"]] <- stats::update.formula(formula, treat ~ .)
  }
  args[["treat"]] <- treat
  args[["covs"]] <- covs
  args[["data"]] <- data
  args[["model"]] <- fit_object
  args[["verbose"]] <- verbose_output
  args[["link"]] <- link



  ################### FITTING THE MODELS #######################################
  if (treat_type == "multinom" || treat_type == "binary" ||
    (treat_type == "ordinal" && method != "polr")) {
    ## --NNET::multinom()-------------------------------------------------------
    if (method == "multinom") {
      infos <- .gps_methods[["multinom"]]
      rlang::check_installed(infos$packages_needed)
      lapply(infos$packages_needed, requireNamespace, quietly = TRUE)

      ## Processing the stuff
      args <- match_add_args(
        arglist = args,
        funlist = lapply(infos$fun.arg.check, function(x) eval(parse(text = x)))
      )

      suppressWarnings(args <- Filter(function(x) !all(is.na(x)), args))

      # fixing single column bug for binary treatments
      if (nunique(args[["data"]][, "treat"]) == 2) {
        former_levels <- levels(args[["data"]][, "treat"])
        levels(args[["data"]][, "treat"]) <- 1:2
        args[["data"]][, "treat"] <- stats::relevel(args[["data"]][, "treat"],
          ref = 1
        )
      }

      ## Fit the multinom
      tryCatch(
        verbosely(
          {
            fit <- do.call(nnet::multinom,
              args = args,
              quote = FALSE
            )
          },
          verbose = verbose_output
        ),
        error = function(e) {
          chk::abort_chk(strwrap(sprintf(
            "There was a problem fitting the multinomial %s regressions
            with `nnet::multinom()`.
            \nError message: (from `nnet::multinom()`) %s",
            infos[["link"]], conditionMessage(e)
          ), prefix = " ", initial = ""), tidy = FALSE)
        }
      )

      # fixing single column bug for binary treatments
      if (nunique(args[["data"]][, "treat"]) == 2) {
        fit_vals <- cbind(1 - fit$fitted.values, fit$fitted.values)
        colnames(fit_vals) <- former_levels
        fit$fitted.values <- fit_vals
      }

      fitted_obj <- fit
    }

    ## --VGLM-------------------------------------------------------------------
    if (method == "vglm") {
      infos <- .gps_methods[["vglm"]]
      rlang::check_installed(infos$packages_needed)
      lapply(infos$packages_needed, requireNamespace, quietly = TRUE)

      ## Procesing additional args
      ## Control
      if (is.null(args[["control"]])) {
        if (link == "multinomial_logit") {
          args[["control"]] <- VGAM::vglm.control()
        } else if (link == "reduced_rank_ml") {
          args[["control"]] <- VGAM::rrvglm.control()
        }
      } else {
        if (link == "multinomial_logit") {
          .chk_cond(
            !args[["control_call"]] ||
              (args[["control_call_char"]] != "VGAM::vglm.control" &&
                args[["control_call_char"]] != "vglm.control"),
            "The argument control has to be a valid function
              call to the function VGAM::vglm.control()"
          )
        } else if (link == "reduced_rank_ml") {
          .chk_cond(
            !args[["control_call"]] ||
              (args[["control_call_char"]] != "VGAM::rrvglm.control" &&
                args[["control_call_char"]] != "rrvglm.control"),
            "The argument control has to be a valid function call
                     to the function VGAM::rrvglm.control()"
          )
        }
      }

      ## family
      if (is.null(args[["family"]])) {
        args[["family"]] <- VGAM::multinomial()
      } else {
        tryCatch(
          {
            VGAM::vglm(args[["formula"]],
              family = args[["family"]],
              data = args[["data"]]
            )
          },
          error = function(e) {
            chk::abort_chk(strwrap("The `family` argument has to be a valid
                                 VGAM family function argument.",
              prefix = " ", initial = ""
            ))
          }
        )
      }

      ## Processing args
      args <- match_add_args(
        arglist = args,
        funlist = lapply(infos$fun.arg.check, function(x) eval(parse(text = x)))
      )

      suppressWarnings(args <- Filter(function(x) !all(is.na(x)), args))

      ## Overwriting args
      ## trace (verbose)
      args[["trace"]] <- verbose_output

      ## method
      if (link == "reduced_rank_ml") args[["method"]] <- "rrvglm.fit"

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
                args = args
              )
            },
            verbose = verbose_output
          ),
          error = function(e) {
            chk::abort_chk(strwrap(sprintf(
              "There was a problem fitting the %s regressions with %s.\n
                               Error message: (from %s) %s",
              link, fun_used, fun_used, conditionMessage(e)
            ), prefix = " ", initial = ""), tidy = FALSE)
          }
        )

        fitted_obj <- fit
      } else {
        probably_a_bug <- TRUE
      }
    }

    ## --brglm2::brmultinom()---------------------------------------------------
    if (method == "brglm2") {
      infos <- .gps_methods[["brglm2"]]
      rlang::check_installed(infos$packages_needed)
      lapply(infos$packages_needed, requireNamespace, quietly = TRUE)

      ## Processin control arg
      if (is.null(args[["control"]])) {
        args[["control"]] <- brglm2::brglmControl()
      } else {
        if (link == "baseline_category_logit") {
          .chk_cond(
            !args[["control_call"]] ||
              (args[["control_call_char"]] != "brglm2::brglmControl" &&
                args[["control_call_char"]] != "brglmControl"),
            "The argument control has to be a valid function call
                    to the function brglm2::brglmControl()"
          )
        }
      }

      ## Processing the arguments
      args <- match_add_args(
        arglist = args,
        funlist = lapply(infos$fun.arg.check, function(x) eval(parse(text = x)))
      )

      suppressWarnings(args <- Filter(function(x) !all(is.na(x)), args))

      ## Fit the brglm2
      if (link %in% infos$link_fun) {
        tryCatch(
          verbosely(
            {
              fit <- do.call(brglm2::brmultinom,
                args = args
              )
            },
            verbose = verbose_output
          ),
          error = function(e) {
            chk::abort_chk(strwrap(sprintf(
              "There was a problem fitting the %s regressions
              with `brglm2::brmultinom()`.\nError message: (from
              `brglm2::brmultinom()`) %s",
              link, conditionMessage(e)
            ), prefix = " ", initial = ""), tidy = FALSE)
          }
        )

        fitted_obj <- fit
      } else {
        probably_a_bug <- TRUE
      }
    }

    ## --mclogit::mclogit()-----------------------------------------------------
    if (method == "mblogit") {
      infos <- .gps_methods[["mblogit"]]
      rlang::check_installed(infos$packages_needed)
      lapply(infos$packages_needed, requireNamespace, quitely = TRUE)
      names.matrix <- unique(args[["treat"]])
      ncol_matrix <- length(names.matrix)

      ## Process the control arg
      if (is.null(args[["control"]])) {
        args[["control"]] <- mclogit::mclogit.control()
      } else {
        if (link == "baseline_category_logit" || "conditional_logit") {
          .chk_cond(
            !args[["control_call"]] ||
              (args[["control_call_char"]] != "mclogit::mclogit.control" &&
                args[["control_call_char"]] != "mclogit.control" &&
                args[["control_call_char"]] != "mclogit::mmclogit.control" &&
                args[["control_call_char"]] != "mmclogit.control"),
            "The argument control has to be a valid function call to the
                  function mclogit::mclogit.control()"
          )
        }
      }

      ## Overwriting args
      if (is.null(args[["estimator"]])) {
        args[["estimator"]] <- "ML"
      }

      ## Processing args
      args <- match_add_args(
        arglist = args,
        funlist = lapply(infos$fun.arg.check, function(x) eval(parse(text = x)))
      )

      suppressWarnings(args <- Filter(function(x) !all(is.na(x)), args))

      if (link %in% infos$link_fun) {
        # Fit model
        tryCatch(
          verbosely(
            {
              fit <- do.call(mclogit::mblogit,
                args = args
              )
            },
            verbose = verbose_output
          ),
          error = function(e) {
            chk::abort_chk(strwrap(sprintf(
              "There was a problem fitting the %s regressions with
              `mclogit::mblogit.\n
              Error message: (from `mclogit::mblogit`) %s",
              link, conditionMessage(e)
            ), prefix = " ", initial = ""), tidy = FALSE)
          }
        )

        fitted_obj <- fit

        ## processing the fitted.values matrix
        fitted.matrix <- matrix(NA,
          nrow = length(fitted_obj$fitted.values) /
            ncol_matrix,
          ncol = ncol_matrix
        )
        colnames(fitted.matrix) <- names.matrix

        for (i in 1:ncol_matrix) {
          sub_vector <- seq(
            from = i,
            to = length(fitted_obj$fitted.values),
            by = ncol_matrix
          )

          fitted.matrix[, i] <- fitted_obj$fitted.values[sub_vector]
        }

        fitted_obj$fitted.values <- fitted.matrix
      } else {
        probably_a_bug <- TRUE
      }
    }

    ## --polr::polr()-----------------------------------------------------------
  } else if (treat_type == "ordinal") {
    if (method == "polr") {
      infos <- .gps_methods[["polr"]]
      rlang::check_installed(infos$packages_needed)

      # map link to method arg of MASS::polr
      which_change <- which(names(args) == "link")
      names(args)[which_change] <- "method"


      ## Processing args
      args <- match_add_args(
        arglist = args,
        funlist = lapply(infos$fun.arg.check, function(x) eval(parse(text = x)))
      )

      suppressWarnings(args <- Filter(function(x) !all(is.na(x)), args))

      if (link %in% infos$link_fun) {
        tryCatch(
          verbosely(
            {
              fit <- do.call(MASS::polr,
                args = args
              )
            },
            verbose = verbose_output
          ),
          error = function(e) {
            chk::abort_chk(strwrap(sprintf(
              "There was a problem fitting the %s regressions with
              `MASS::polr()`.\nError message: (from `MASS::polr()`) %s",
              link, conditionMessage(e)
            ), prefix = " ", initial = ""), tidy = FALSE)
          }
        )
      }
      fitted_obj <- fit
    }
  } else {
    probably_a_bug <- TRUE
  }

  ## Last check of output object
  .chk_cond(
    probably_a_bug,
    "The function `estimate_gps()` was not able to estimate the
            propensity scores.  It's probably a bug. Please let the author
            know."
  )

  return(fitted_obj)
}
