## --extracting varaibles from a formula----------------------------------------
## --based on: https://github.com/ngreifer/WeightIt/blob/master/R/utils.R-------

#' @srrstats {G1.4a} All helpers have roxygen2-style docstrings
#' Internal helper to extract variables from a formula
#'
#' Parses a model formula and associated data to extract the treatment variable,
#' reported covariates, and (optionally) a model matrix including interactions.
#' Handles formulas with data.frame terms, checks that all variables exist in
#' `data` or the parent environment, and prevents the treatment from appearing
#' on the right-hand side. Adapted from utilities in the WeightIt package.
#' @noRd
.get_formula_vars <- function(formula, data = NULL, ...) {
  parent_env <- environment(formula)
  eval_model_matrx <- !(any(c("|", "||") %in% all.names(formula)))

  ## Check data; if not exists the look for the data in the parent env----------
  if (!is.null(data)) {
    .check_df(data)
  } else {
    data <- environment(formula)
  }

  ## Check formula--------------------------------------------------------------
  if (missing(formula) || !rlang::is_formula(formula)) {
    chk::abort_chk(strwrap("The argument formula has to be a valid R formula",
      prefix = " ", initial = ""
    ))
  }

  ## Extract stats::terms from the formula
  tryCatch(vars <- stats::terms(formula, data = data),
    error = function(e) {
      chk::abort_chk(conditionMessage(e), tidy = TRUE)
    }
  )

  ## extract treatment from the data or env
  if (rlang::is_formula(vars, lhs = TRUE)) {
    treat <- attr(vars, "variables")[[2]]
    treat_name <- deparse(treat)

    tryCatch(treat_data <- eval(treat, data, parent_env),
      error = function(e) {
        chk::abort_chk(conditionMessage(e), tidy = TRUE)
      }
    )
  } else {
    chk::abort_chk(strwrap("The argument formula has to be a valid R formula",
      prefix = " ", initial = ""
    ))
  }

  ## handle covariates - right-hand-side variables (RHS)------------------------
  vars_covs <- stats::delete.response(vars)
  rhs_vars <- attr(vars_covs, "variables")[-1]
  rhs_vars_char <- vapply(rhs_vars, deparse, character(1L))

  covs <- list()
  lapply(seq_along(rhs_vars), function(i) {
    tryCatch(eval(rhs_vars[[i]], data, parent_env),
      error = function(e) {
        chk::abort_chk(strwrap("All variables in the `formula` must be columns
                               in the `data` or objects in the global
                               environment.", prefix = " ", initial = ""))
      }
    )
  })

  ## dealing with interactions--------------------------------------------------
  rhs_labels <- attr(vars_covs, "term.labels")
  rhs_labels_list <- stats::setNames(as.list(rhs_labels), rhs_labels)
  rhs_order <- attr(vars_covs, "order")

  ## additional dfs in the formula
  rhs_if_df <- stats::setNames(vapply(rhs_vars, function(v) {
    length(dim(try(eval(v, data, parent_env)))) == 2L
  }, logical(1L)), rhs_vars_char)

  if (any(rhs_if_df)) {
    if (any(rhs_vars_char[rhs_if_df] %in%
      unlist(lapply(
        rhs_labels[rhs_order > 1],
        function(x) strsplit(x, ":", fixed = TRUE)
      )))) {
      chk::abort_chk("Interactions with data.frames are not allowed.")
    }

    addl_dfs <- stats::setNames(
      lapply(which(rhs_if_df), function(i) {
        df <- eval(rhs_vars[[i]], data, parent_env)
        if (inherits(df, "rms")) {
          class(df) <- "matrix"
          df <- stats::setNames(
            as.data.frame(as.matrix(df)),
            attr(df, "colnames")
          )
        } else {
          colnames(df) <- paste(rhs_vars_char[i], colnames(df), sep = "_")
        }
        df <- as.data.frame(df)
      }),
      rhs_vars_char[rhs_if_df]
    )

    for (i in rhs_labels[rhs_labels %in% rhs_vars_char[rhs_if_df]]) {
      ind <- which(rhs_labels == i)
      rhs_labels <- append(rhs_labels[-ind],
        values = names(addl_dfs[[i]]),
        after = ind - 1
      )
      rhs_labels_list[[i]] <- names(addl_dfs[[i]])
    }

    if (!is.null(data)) {
      data <- do.call("cbind", unname(c(addl_dfs, list(data))))
    } else {
      data <- do.call("cbind", unname(addl_dfs))
    }
  }

  ## dealing with no stats::terms-----------------------------------------------
  if (is.null(rhs_labels)) {
    new_form <- stats::as.formula("~0")
    vars_covs <- stats::terms(new_form)
    covs <- data.frame(Intercept = rep.int(1, if (is.null(treat)) {
      1L
    } else {
      length(treat)
    }))[, -1, drop = FALSE]
  } else {
    new_form_char <- sprintf("~ %s", paste(vapply(
      names(rhs_labels_list), function(x) {
        if (x %in% rhs_vars_char[rhs_if_df]) {
          paste0("`", rhs_labels_list[[x]], "`", collapse = " + ")
        } else {
          rhs_labels_list[[x]]
        }
      },
      character(1L)
    ), collapse = " + "))

    new_form <- stats::as.formula(new_form_char)
    vars_covs <- stats::terms(stats::update(new_form, ~ . - 1))

    # Get model.frame
    mf_covs <- quote(stats::model.frame(vars_covs, data,
      drop.unused.levels = TRUE,
      na.action = "na.pass"
    ))

    tryCatch(
      {
        covs <- eval(mf_covs)
      },
      error = function(e) {
        chk::abort_chk(conditionMessage(e), tidy = TRUE)
      }
    )

    if (!is.null(treat_name) && treat_name %in% names(covs)) {
      chk::abort_chk(strwrap("The treatment variable cannot appear on the right
      side of the formula", prefix = " ", initial = ""))
    }
  }

  if (eval_model_matrx) {
    original_covs_levels <- .make_list(names(covs))

    for (i in names(covs)) {
      if (is.character(covs[[i]])) {
        covs[[i]] <- factor(covs[[i]])
      } else if (!is.factor(covs[[i]])) {
        next
      }

      if (length(unique(covs[[i]])) == 1L) {
        covs[[i]] <- 1
      } else {
        original_covs_levels[[i]] <- levels(covs[[i]])
        levels(covs[[i]]) <- paste0("", original_covs_levels[[i]])
      }
    }

    # Get full model matrix with interactions too
    covs_matrix <- stats::model.matrix(vars_covs,
      data = covs,
      contrasts.arg = lapply(Filter(is.factor, covs),
        stats::contrasts,
        contrasts = FALSE
      )
    )

    for (i in names(covs)[vapply(covs, is.factor, logical(1L))]) {
      levels(covs[[i]]) <- original_covs_levels[[i]]
    }
  } else {
    covs_matrix <- NULL
  }

  # Defining the list to return-------------------------------------------------
  ret <- list(
    treat = treat_data,
    treat_name = treat_name,
    reported_covs = covs,
    model_covs = covs_matrix
  )
  return(ret)
}

#' Internal helper to process model formulas for GPS estimation
#'
#' Validates the input formula, ensuring it has both treatment and predictor
#' variables, and uses `.get_formula_vars()` to extract the treatment and
#' covariate model matrix. Performs basic consistency checks (length, missing
#' values, and number of treatment levels) and warns when the treatment has
#' many levels before returning the parsed components.
#' @noRd
.process_formula <- function(formula, data) {
  # defining empty list for arg storage
  args <- list()

  # formula
  if (missing(formula)) {
    chk::abort_chk(strwrap("The argument `formula` is missing with no default",
      prefix = " ", initial = ""
    ))
  }

  if (!rlang::is_formula(formula, lhs = TRUE)) {
    chk::abort_chk(strwrap("The argument `formula` has to be a valid R formula
    with treatment and predictor variables", prefix = " ", initial = ""))
  }

  data_list <- .get_formula_vars(formula, data)

  args["treat"] <- list(data_list[["treat"]])
  args["covs"] <- list(data_list[["model_covs"]])

  if (is.null(args["treat"])) {
    chk::abort_chk(strwrap("No treatment variable was specified",
      prefix = " ", initial = ""
    ))
  }

  if (is.null(args["covs"])) {
    chk::abort_chk(strwrap("No predictors were specified",
      prefix = " ", initial = ""
    ))
  }

  if (length(args[["treat"]]) != nrow(args[["covs"]])) {
    chk::abort_chk(strwrap("The treatment variable and predictors ought to have
    the same number of samples", prefix = " ", initial = ""))
  }

  if (anyNA(args[["treat"]])) {
    chk::abort_chk(strwrap("The `treatment` variable can not have any NA's",
      prefix = " ", initial = ""
    ))
  }

  n_levels <- nunique(args[["treat"]])
  if (n_levels > 10) {
    .vm_warn(strwrap("The `treatment` variable has more than 10 unique levels.
    Consider dropping the number of groups, as the vector matching algorithm may
             not perform well", prefix = " ", initial = ""))
  }

  # return list with output vars
  return(data_list)
}

#' Internal helper to construct lists
#'
#' Creates a list of the appropriate length, optionally named from `n`.
#' @noRd
.make_list <- function(n) {
  if (length(n) == 1L && is.numeric(n)) {
    vector("list", as.integer(n))
  } else if (length(n) > 0L && is.atomic(n)) {
    stats::setNames(vector("list", length(n)), as.character(n))
  } else {
    chk::abort_chk(strwrap("'n' must be an integer(ish) scalar or an atomic
                           variable.", prefix = " ", initial = ""))
  }
}

#' Internal helper to count unique values
#'
#' Returns the number of unique values in a vector, with special handling for
#' factors and optional removal of missing values.
#' @noRd
nunique <- function(x, na_rm = TRUE) {
  if (is.null(x)) {
    return(0)
  }
  if (is.factor(x)) {
    return(nlevels(x))
  }
  if (na_rm && anyNA(x)) x <- na_rem(x)
  length(unique(x))
}

#' Internal helper to remove missing values
#'
#' A simple, fast `na.omit()` alternative for vectors.
#' @noRd
na_rem <- function(x) {
  # A faster na.omit for vectors
  x[!is.na(x)]
}

#' Internal helper to format word lists
#'
#' Builds human-readable lists for error and message text (e.g. "a and b" or
#' "a, b, and c"), with optional quoting and verb agreement.
#' @noRd
word_list <- function(word_list = NULL, and_or = "and", is_are = FALSE,
                      quotes = FALSE) {
  # When given a vector of strings, creates a string of the form "a and b"
  # or "a, b, and c"
  # If is_are, adds "is" or "are" appropriately

  word_list <- setdiff(word_list, c(NA_character_, ""))

  if (is.null(word_list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word_list <- add_quotes(word_list, quotes)

  len_wl <- length(word_list)

  if (len_wl == 1L) {
    out <- word_list
    if (is_are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is.null(and_or) || isFALSE(and_or)) {
    out <- paste(word_list, collapse = ", ")
  } else {
    and_or <- match.arg(and_or, c("and", "or"))

    if (len_wl == 2L) {
      out <- sprintf(
        "%s %s %s",
        word_list[1L],
        and_or,
        word_list[2L]
      )
    } else {
      out <- sprintf(
        "%s, %s %s",
        paste(word_list[-len_wl], collapse = ", "),
        and_or,
        word_list[len_wl]
      )
    }
  }

  if (is_are) out <- sprintf("%s are", out)

  attr(out, "plural") <- TRUE

  out
}

#' Internal helper to add quotes around strings
#'
#' Wraps character vectors in single or double quotes, or leaves them unchanged,
#' depending on the `quotes` argument.
#' @noRd
add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes)) {
    quotes <- '"'
  }

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1) {
      sprintf("'%s'", x)
    } else {
      sprintf('"%s"', x)
    }
  }

  x
}

#' Internal helper to process stratification variable
#'
#' Validates and processes the `by` argument used for stratifying matching or
#' estimation, extracting the grouping column from `data`, checking for missing
#' values, and ensuring that each stratum contains all treatment levels. Returns
#' a one-column data frame with an attached `by.factor` attribute encoding the
#' strata.
#' @noRd
.process_by <- function(by, data, treat) {
  ## Process by
  error_by <- FALSE
  n <- length(treat)

  if (missing(by)) {
    error_by <- TRUE
  } else if (is.null(by)) {
    return(NULL)
  } else if (chk::vld_string(by) && by %in% colnames(data)) {
    by.data <- data[[by]]
    by.name <- by
  } else {
    error_by <- TRUE
  }

  .chk_cond(
    error_by,
    "The argument `by` must be a single string with the name of
                   the column to stratify by, or a one sided formula with one
                   stratifying variable on the right-hand side"
  )

  .chk_cond(
    anyNA(by.data),
    "The argument `by` cannot contain any NA's"
  )

  by_comps <- data.frame(by.data)

  names(by_comps) <- {
    if (!is.null(colnames(by.data))) {
      colnames(by.data)
    } else {
      by.name
    }
  }

  by.factor <- {
    if (is.null(by.data)) {
      factor(rep.int(1L, n), levels = 1L)
    } else {
      factor(by_comps[[1]],
        levels = sort(unique(by_comps[[1]])),
        labels = paste0(names(by_comps), "=", sort(unique(by_comps[[1]])))
      )
    }
  }

  if (any(vapply(
    levels(by.factor),
    function(x) nunique(treat) != nunique(treat[by.factor == x]),
    logical(1L)
  ))) {
    chk::abort_chk(strwrap("Not all group formed by the column provided to the
    argument `by` contain all treatment levels.", prefix = " ", initial = ""))
  }

  attr(by_comps, "by.factor") <- by.factor

  return(by_comps)
}

## --scaling the data-----------------------------------------------------------
#' Internal helper to test for binary variables
#'
#' Checks whether a vector is binary, either restricted to 0/1 or allowing any
#' two distinct values.
#' @noRd
is_binary <- function(x, null_one = TRUE) {
  if (null_one == TRUE) {
    return(all(stats::na.omit(x) %in% 0:1))
  } else {
    return(length(unique(x)) == 2L)
  }
}

#' Internal helper to test if all values are identical
#'
#' Checks whether all (non-missing) elements of a vector are the same.
#' @noRd
all_the_same <- function(x, na_rm = TRUE) {
  if (anyNA(x)) {
    x <- x[!is.na(x)]
    if (!na_rm) {
      return(is.null(x))
    }
  }

  if (is.numeric(x)) {
    (max(x) - min(x)) == 0L
  } else {
    all(x == x[1])
  }
}

#' Internal helper to rescale variables to (0, 1)
#'
#' Leaves binary and constant variables unchanged, rescales numeric variables
#' to (0, 1), and coerces two-level variables to factors.
#' @noRd
scale_0_to_1 <- function(x) {
  if (is_binary(x)) {
    return(x)
  } else if (is_binary(x, null_one = FALSE)) {
    if (is.character(x)) {
      as.factor(x)
    } else if (is.logical(x)) {
      x
    } else {
      is_first <- {
        x == unique(x)[1]
      }
      x <- ifelse(is_first, 0, 1)
      return(as.factor(x))
    }
  } else if (is.factor(x) || all_the_same(x)) {
    return(x)
  } else if (is.numeric(x)) {
    if (min(x) >= 0 && max(x) <= 1) {
      return(x)
    } else {
      x_scaled <- (x - min(x)) / (max(x) - min(x))
      return(x_scaled)
    }
  } else {
    chk::abort_chk(strwrap("Invalid type. Cannot convert x to 0-1 range"),
      prefix = " ", initial = ""
    )
    return(as.factor(x))
  }
}

`%nin%` <- function(x, inx) {
  !(x %in% inx)
}

#' Internal helper to align argument lists with function formals
#'
#' Matches names in an argument list to the formals of one or more functions,
#' dropping unused entries and filling in defaults where needed.
#' @noRd
match_add_args <- function(arglist, funlist) {
  if (is.list(arglist) && length(arglist) == sum(names(arglist) != "",
    na.rm = TRUE
  )) {
    argnames <- names(arglist)

    formlist <- {
      if (length(funlist) == 1L && is.function(funlist)) {
        formals(funlist)
      } else if (is.list(funlist)) {
        unlist(lapply(funlist, function(x) {
          forms <- formals(x)
          forms[unlist(lapply(forms, is.null))] <- NA
          forms
        }))
      }
    }

    ## Check for doubles and resolve to default if present
    formnames <- names(formlist)

    if (length(funlist) != 1 && !is.function(funlist)) {
      duplicates <- formnames[duplicated(formnames)]

      remove_forms <- numeric()
      for (i in seq_along(duplicates)) {
        cur_dups <- which(formnames %in% duplicates[i])

        is_empty <- lapply(formlist[cur_dups], function(x) {
          is.symbol(x) && as.character(x) == ""
        })

        if (all(unlist(is_empty)) || all(!unlist(is_empty))) {
          remove_forms <- c(remove_forms, cur_dups[-1])
        } else {
          not_empty <- cur_dups[!unlist(is_empty)]
          remove_forms <- c(remove_forms, not_empty[-1])
        }
      }

      formlist <- formlist[-remove_forms]
      formnames <- names(formlist)
    }

    ## match the names
    matched_names <- vapply(argnames, function(x) {
      match(x, formnames)
    }, numeric(1))

    ## Set names of present arguments
    choose_names <- names(matched_names[!is.na(matched_names)])
    arglist <- arglist[names(arglist) %in% choose_names]

    ## Add default values of non defined args
    which_undefined <- which(formnames %nin% argnames)

    is_empty2 <- lapply(formlist[which_undefined], function(x) {
      is.symbol(x) && as.character(x) == ""
    })

    if (length(is_empty2) != 0) {
      add_formals <- which_undefined[!unlist(is_empty2)]
      arglist <- append(arglist, formlist[add_formals])
    } else {
      arglist <- append(arglist, formlist)
    }

    return(arglist)
  }
}

#' Internal helper to classify treatment type
#'
#' Converts the treatment variable to a factor, checks variability, and assigns
#' it a type label (`"binary"`, `"ordinal"`, or `"multinom"`) via the
#' `treat_type` attribute.
#' @noRd
.assign_treatment_type <- function(treat, ordinal_treat) {
  chk::chk_vector(treat)

  treat <- as.factor(treat)

  if (all_the_same(treat)) {
    chk::abort_chk(strwrap("There is no variability in the treatment variable.
    All datapoints are the same.", prefix = " ", initial = ""))
  } else if (is_binary(treat, null_one = FALSE)) {
    treat_type <- "binary"
  } else if (is.factor(treat) && is.ordered(treat)) {
    treat_type <- "ordinal"
  } else {
    treat_type <- "multinom"
  }

  attr(treat, "treat_type") <- treat_type
  treat
}

#' Internal helper to retrieve treatment type
#'
#' Returns the `treat_type` attribute from a treatment vector.
#' @noRd
.get_treat_type <- function(treat) {
  attr(treat, "treat_type")
}

verbosely <- function(expr, verbose = TRUE) {
  if (verbose) {
    return(expr)
  }

  utils::capture.output({
    out <- invisible(expr)
  })

  out
}

## --get rid of cli package note------------------------------------------------
ignore_unused_imports <- function() {
  c(
    cli::cli_warn,
    Matching::Match,
    optmatch::match_on
  )
}

#' Internal helper to process reference level
#'
#' Validates and sets the reference level of the treatment variable, optionally
#' warning when `reference` is redundant given `ordinal_treat`. Returns the
#' relevelled factor and the reference used.
#' @noRd
.process_ref <- function(data_vec,
                         reference = NULL,
                         ordinal_treat = NULL) {
  # reference
  levels_treat <- as.character(unique(data_vec))

  if (!is.null(ordinal_treat) && !is.null(reference)) {
    .vm_warn(strwrap("There is no need to specify `reference` if `ordinal_treat`
                     was provided. Ignoring the `reference` argument",
      prefix = " ", initial = ""
    ))
  } else {
    if (is.null(reference)) {
      reference <- levels_treat[1]

      if (!is.ordered(data_vec)) {
        data_vec <- stats::relevel(data_vec, ref = reference)
      }
    } else if (!(is.character(reference) && length(reference) == 1L &&
      !anyNA(reference))) {
      chk::abort_chk(strwrap("The argument `reference` must be a single string
                             of length 1", prefix = " ", initial = ""))
    } else if (reference %nin% levels_treat) {
      chk::abort_chk(strwrap("The argument `reference` is not in the unique
      levels of the treatment variable", prefix = " ", initial = ""))
    } else {
      data_vec <- factor(data_vec, ordered = FALSE)
      data_vec <- stats::relevel(data_vec, ref = reference)
    }
  }

  ## assembling output
  ref_out <- list(
    data.relevel = data_vec,
    reference = reference
  )

  return(ref_out)
}

#' Internal helper for logit transformation
#'
#' Applies the logit transform \eqn{\log(x / (1 - x))} to probabilities in
#' (0, 1).
#' @noRd
logit <- function(x) {
  log(x / (1 - x))
}

#' Internal helper to vectorize scalar arguments
#'
#' Recycles a length-1 argument to the requested length, otherwise returns it
#' unchanged.
#' @noRd
.vectorize <- function(arg, times) {
  if (length(arg) == 1) {
    rep(arg, times)
  } else {
    arg
  }
}

#' Internal helper to summarize balance quality tables
#'
#' Aggregates balance statistics across time, reshapes them, and classifies
#' coefficients as balanced or not for use in `balqual()` output.
#' @noRd
create_balqual_output <- function(coeflist,
                                  coefnames,
                                  operation = c("+", "max"),
                                  round = 2,
                                  which_coefs = "smd",
                                  cutoffs = c(0.1, 0.1, 2)) {
  # summarize the coeflist
  if (operation == "+") {
    quality_table <- Reduce("+", coeflist) / length(coeflist)
  } else if (operation == "max") {
    quality_table <- do.call(pmax, coeflist)
  }

  # round the output
  quality_table <- round(quality_table, round)

  # cbind with coefficient names
  quality_table <- cbind(coefnames, quality_table)

  # stats::reshape variables to long
  quality_table <- stats::reshape(
    quality_table,
    direction = "long",
    varying = list(names(quality_table)[3:ncol(quality_table)]),
    v.names = "value",
    idvar = c("coef_name", "time"),
    timevar = "variable",
    times = names(quality_table)[3:ncol(quality_table)]
  )

  # stats::reshape the time to wide
  quality_table <- stats::reshape(
    quality_table,
    idvar = c("coef_name", "variable"),
    timevar = "time",
    direction = "wide"
  )

  # get rid of the rownames and arrange the output
  rownames(quality_table) <- NULL
  quality_table <- quality_table[, c(
    "variable", "coef_name", "value.before",
    "value.after"
  )]


  # calculate the quality
  quality_table$quality <- rep(NA, nrow(quality_table))

  filter_list <- lapply(list("smd", "r", "var_ratio"), function(x) {
    quality_table$coef_name == x
  })

  bal_values <- mapply(
    function(x, y) {
      ifelse(
        quality_table[x, "value.after"] < y,
        "Balanced",
        "Not Balanced"
      )
    },
    filter_list,
    cutoffs,
    SIMPLIFY = FALSE
  )

  for (i in seq_len(length(filter_list))) {
    quality_table$quality[filter_list[[i]]] <- bal_values[[i]]
  }

  # changing colnames
  colnames(quality_table) <- c(
    "Variable", "Coefficient", "Before", "After",
    "Reduction"
  )

  filter_rows <- quality_table$Coefficient %in% which_coefs

  # change coef recoding
  quality_table$Coefficient <- factor(quality_table$Coefficient,
    levels = c("smd", "r", "var_ratio"),
    labels = c("SMD", "r", "Var")
  )

  # filter output
  quality_table[filter_rows, ]
}

#' Internal helper to order GPS vectors
#'
#' Produces indices ordering a GPS vector in descending, ascending, original, or
#' random order for use in matching.
#' @noRd
.ordering_func <- function(gps_vec, order = "desc", ...) {
  if (order == "desc") {
    order(gps_vec, decreasing = TRUE, ...)
  } else if (order == "asc") {
    order(gps_vec, decreasing = FALSE, ...)
  } else if (order == "original") {
    seq_len(length(gps_vec))
  } else if (order == "random") {
    withr::with_preserve_seed(sample(gps_vec))
  }
}

# declare global variables to hanlde data-mask notes
.onLoad <- function(libname, pkgname) {
  # suppress spurious no-visible-binding notes for data-mask vars
  utils::globalVariables(c("gps_method", "min_controls", "max_controls", "i"))
}
