##--extracting varaibles from a formula-----------------------------------------
##--based on: https://github.com/ngreifer/WeightIt/blob/master/R/utils.R--------

.get_formula_vars <- function(formula, data = NULL, ...) {
  args <- list(...)
  parent_env <- environment(formula)
  eval.model.matrx <- !(any(c("|", "||") %in% all.names(formula)))

  ## Check data; if not exists the look for the data in the parent env----------
  if(!is.null(data)) {
    .check_df(data)
  } else {
    data <- environment(formula)
  }

  ## Check formula--------------------------------------------------------------
  if(missing(formula) || !rlang::is_formula(formula)) {
    chk::abort_chk('The argument formula has to be a valid R formula')
  }

  ## Extract terms from the formula
  tryCatch(vars <- stats::terms(formula, data = data),
           error = function(e) {
             chk::abort_chk(conditionMessage(e), tidy = TRUE)
           })

  ## extract treatment from the data or env
  if(rlang::is_formula(vars, lhs = TRUE)) {
    treat <- attr(vars, 'variables')[[2]]
    treat_name <- deparse(treat)

    tryCatch(treat_data <- eval(treat, data, parent_env),
             error = function(e) {
               chk::abort_chk(conditionMessage(e), tidy = TRUE)
             })
  } else {
    chk::abort_chk('The argument formula has to be a valid R formula')
  }

  ## handle covariates - right-hand-side variables (RHS)------------------------
  vars_covs <- delete.response(vars)
  rhs_vars <- attr(vars_covs, 'variables')[-1]
  rhs_vars_char <- vapply(rhs_vars, deparse, character(1L))

  covs <- list()
  lapply(seq_along(rhs_vars), function(i) {
    tryCatch(test <- eval(rhs_vars[[i]], data, parent_env), #covs[[i]] <<-
             error = function(e) {
               chk::abort_chk('All variables in the `formula` must be columns
                              in the `data` or objects in the global environment.')
             })
  })

  #covs <- setNames(covs, rhs_vars_char)

  ## dealing with interactions--------------------------------------------------
  rhs_labels <- attr(vars_covs, 'term.labels')
  rhs_labels_list <- setNames(as.list(rhs_labels), rhs_labels)
  rhs_order <- attr(vars_covs, 'order')

  ## additional dfs in the formula
  rhs_if_df <- setNames(vapply(rhs_vars, function(v) {
    length(dim(try(eval(v, data, parent_env)))) == 2L
  }, logical(1L)), rhs_vars_char)

  if(any(rhs_if_df)) {
    if(any(rhs_vars_char[rhs_if_df] %in%
           unlist(lapply(rhs_labels[rhs_order > 1],
                         function(x) strsplit(x, ':', fixed = TRUE))))) {
      chk::abort_chk('Interactions with data.frames are not allowed.')
    }

    addl_dfs <- setNames(lapply(which(rhs_if_df), function(i) {
      df <- eval(rhs_vars[[i]], data, parent_env)
      if (inherits(df, 'rms')) {
        class(df) <- 'matrix'
        df <- setNames(as.data.frame(as.matrix(df)), attr(df, 'colnames'))
      } else {
        colnames(df) <- paste(rhs_vars_char[i], colnames(df), sep = '_')
      }
      df <- as.data.frame(df)
    }),
    rhs_vars_char[rhs_if_df])

    for(i in rhs_labels[rhs_labels %in% rhs_vars_char[rhs_if_df]]) {
      ind <- which(rhs_labels == i)
      rhs_labels <- append(rhs_labels[-ind],
                           values = names(addl_dfs[[i]]),
                           after = ind - 1)
      rhs_labels_list[[i]] <- names(addl_dfs[[i]])
    }

    if(!is.null(data)) data <- do.call('cbind', unname(c(addl_dfs, list(data))))
    else data <- do.call('cbind', unname(addl_dfs))
  }

  ## dealing with no terms------------------------------------------------------
  if(is.null(rhs_labels)) {
    new_form <- as.formula('~0')
    vars_covs <- terms(new_form)
    covs <- data.frame(Intercept = rep.int(1, if (is.null(treat)) 1L else
      length(treat)))[, -1, drop = FALSE]

  } else {
    new_form_char <- sprintf("~ %s", paste(vapply(names(rhs_labels_list), function(x) {
      if (x %in% rhs_vars_char[rhs_if_df]) paste0("`", rhs_labels_list[[x]], "`", collapse = " + ")
      else rhs_labels_list[[x]]},
      character(1L)), collapse = " + "))

    new_form <- as.formula(new_form_char)
    vars_covs <- terms(update(new_form,  ~ . - 1))

    # Get model.frame
    mf_covs <- quote(stats::model.frame(vars_covs, data,
                                        drop.unused.levels = TRUE,
                                        na.action = 'na.pass'))

    tryCatch({covs <- eval(mf_covs)},
             error = function(e) {
               chk::abort_chk(conditionMessage(e), tidy = TRUE)
             })

    if(!is.null(treat_name) && treat_name %in% names(covs)) {
      chk::abort_chk('The treatment variable cannot appear on the right side
                     of the formula')
    }
  }

  if(eval.model.matrx) {
    original_covs_levels <- .make_list(names(covs))

    for (i in names(covs)) {
      if (is.character(covs[[i]]))
        covs[[i]] <- factor(covs[[i]])
      else if (!is.factor(covs[[i]]))
        next

      if (length(unique(covs[[i]])) == 1L) {
        covs[[i]] <- 1
      }
      else {
        original_covs_levels[[i]] <- levels(covs[[i]])
        levels(covs[[i]]) <- paste0("", original_covs_levels[[i]])
      }
    }

    #Get full model matrix with interactions too
    covs_matrix <- model.matrix(vars_covs, data = covs,
                                contrasts.arg = lapply(Filter(is.factor, covs),
                                                       contrasts, contrasts = FALSE))

    for (i in names(covs)[vapply(covs, is.factor, logical(1L))]) {
      levels(covs[[i]]) <- original_covs_levels[[i]]
    }

  } else {
    covs_matrix <- NULL
  }

  # Defining the list to return-------------------------------------------------
  ret <- list(treat = treat_data,
              treat_name = treat_name,
              reported_covs = covs,
              model_covs = covs_matrix)
  return(ret)
}


#R Processing-------------------------------------------------------------------
.make_list <- function(n) {
  if (length(n) == 1L && is.numeric(n)) {
    vector("list", as.integer(n))
  }
  else if (length(n) > 0L && is.atomic(n)) {
    setNames(vector("list", length(n)), as.character(n))
  }
  else stop("'n' must be an integer(ish) scalar or an atomic variable.")
}

#Uniqueness---------------------------------------------------------------------
nunique <- function(x, na.rm = TRUE) {
  if (is.null(x)) return(0)
  if (is.factor(x)) return(nlevels(x))
  if (na.rm && anyNA(x)) x <- na.rem(x)
  length(unique(x))
}

na.rem <- function(x) {
  #A faster na.omit for vectors
  x[!is.na(x)]
}

##--wordlists for error generation----------------------------------------------
word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately

  word.list <- setdiff(word.list, c(NA_character_, ""))

  if (is.null(word.list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word.list <- add_quotes(word.list, quotes)

  L <- length(word.list)

  if (L == 1L) {
    out <- word.list
    if (is.are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is.null(and.or) || isFALSE(and.or)) {
    out <- paste(word.list, collapse = ", ")
  }
  else {
    and.or <- match.arg(and.or, c("and", "or"))

    if (L == 2L) {
      out <- sprintf("%s %s %s",
                     word.list[1L],
                     and.or,
                     word.list[2L])
    }
    else {
      out <- sprintf("%s, %s %s",
                     paste(word.list[-L], collapse = ", "),
                     and.or,
                     word.list[L])
    }
  }

  if (is.are) out <- sprintf("%s are", out)

  attr(out, "plural") <- TRUE

  out
}

add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes))
    return(x)

  if (isTRUE(quotes))
    quotes <- '"'

  if (chk::vld_string(quotes))
    return(paste0(quotes, x, quotes))

  if (!chk::vld_count(quotes) || quotes > 2) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L)
    return(x)

  x <- {
    if (quotes == 1) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}
expand_grid_string <- function(..., collapse = "") {
  do.call("paste", c(expand.grid(...), sep = collapse))
}

##--processing `missing` argument-----------------------------------------------
.process_missing <- function(missing, method) {
  if(is.null(method)) {
    return('')
  }

  allowable.missings <- .gps_methods[[method]]$missing

  chk::chk_null_or(missing, chk::chk_string)

  if(is.null(missing) || !(missing %in% allowable.missings)) {

    chk::abort_chk(sprintf('Only %s allowed for the argument `missing` with
                           the method: %s', word_list(allowable.missings,
                                                      quotes = 2, is.are = TRUE),
                           method))
    return(allowable.missings[1])
  }
}

##--processing `by` argument----------------------------------------------------
.process_by <- function(by, data, treat) {

  ##Process by
  error.by <- FALSE
  n <- length(treat)

  if(missing(by)) {
    error.by <- TRUE
  } else if(is.null(by)) {
    by <- NULL
  } else if(chk::vld_string(by) && by %in% colnames(data)) {
    by <- data[[by]]
  } else if(length(dim(by)) == 2L && nrow(by) == n) {
    by <- drop(by[, 1])
  } else if(rlang::is_formula(by, lhs = FALSE)) {
    covs <- .get_formula_vars(formula, data)
    by <- covs[['reported.covs']]

    if(ncol(by) != 1L) {
      chk::abort_chk('The formula used in the `by` argument is only allowed to
      have variable on the right side ')
    }
  } else {
    error.by <- TRUE
  }

  if(error.by) {
    chk::abort_chk()
  }

}
