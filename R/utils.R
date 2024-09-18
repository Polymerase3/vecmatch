##--extracting varaibles from a formula-----------------------------------------
##--based on: https://github.com/ngreifer/WeightIt/blob/master/R/utils.R--------

.get_formula_vars <- function(formula, data = NULL, ...) {
  args <- list(...)
  parent_env <- environment(formula)

  ## Check data; if not exists the look for the data in the parent env
  if(!is.null(data)) {
    .check_df(data)
  } else {
    data <- environment(formula)
  }

  ## Check formula
  if(missing(formula) || !rlang::is_formula(formula)) {
    chk::abort_chk('The argument formula has to be a valid R formula')
  }

  ## Extract terms from the formula
  tryCatch(vars <- stats::terms(formula, data = data),
           error = function(e) {
             chk::abort_chk(conditionMessage(e), tidy = TRUE)
           })

  ## extract treatment from the data or env


}
