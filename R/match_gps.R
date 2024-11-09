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
match_gps <- function(formula,
                      data = NULL,
                      caliper = NULL,
                      reference = NULL,
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
  if (!is.null(data)) .check_df(data)

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


}
