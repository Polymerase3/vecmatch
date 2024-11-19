#' Title
#'
#' @param formula
#' @param data
#' @param caliper
#' @param reference
#' @param combos
#' @param kmeans.args
#' @param ratio
#' @param replace
#' @param ...
#'
#' @return
#'
#' @examples
#' @export
match_gps <- function(treatment = NULL,
                      data = NULL,
                      caliper = NULL,
                      reference = NULL,
                      by = NULL,
                      combos = NULL,
                      kmeans.args = NULL,
                      ratio = 1,
                      replace = FALSE,
                      ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  args <- list(...)

  # data
  if (!is.null(data)) {
    .check_df(data)
    if('gps' %nin% class(data)) {
      'The argument `data` has to be of class `gps`.'
    }
  } else {
    chk::abort_chk('The argument `data` is missing with no default.')
  }

  # treatment
  symlist <- list(
    treatment = substitute(treatment)
  )

  symlist <- .conv_nam(symlist)

  # check data dimensions
  treat_levels <- unique(data[, symlist[['treatment']]])
  print(treat_levels)





  # formula
  data.list <- .process_formula(formula, data)

  # args assignment to list used in calculations
  args["treat"] <- list(data.list[["treat"]])
  args["covs"] <- list(data.list[["model_covs"]])

  # reference
  ref.list <- .process_ref(args[['treat']],
                           ordinal.treat = NULL,
                           reference = reference)

  args[['treat']] <- ref.list[['data.relevel']]
  reference <- ref.list[['reference']]

  # process caliper
  if(is.null(caliper)) {

  } else {
    chk::chk_numeric(caliper)
    chk::chk_range(caliper, range = c(0, 1))
  }

  ######################## MATCHING ############################################


}
