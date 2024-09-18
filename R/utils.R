##--extracting varaibles from a formula-----------------------------------------
##--based on: https://github.com/ngreifer/WeightIt/blob/master/R/utils.R--------

.get_formula_vars <- function(formula, data = NULL, ...) {
  args <- list(...)
  parent_env <- environment(formula)

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
    tryCatch(covs[[i]] <<- eval(rhs_vars[[i]], data, parent_env),
             error = function(e) {
               chk::abort_chk('All variables in the `formula` must be columns
                              in the `data` or objects in the global environment.')
             })
  })

  covs <- setNames(covs, rhs_vars_char)

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
  }


  # Defining the list to return-------------------------------------------------
  ret <- list(treat = treat_data,
              treat_name = treat_name,
              covs = covs)
  #return(ret)
}
