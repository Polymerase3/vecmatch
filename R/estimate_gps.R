
#' Title
#'
#' @return
#' @export
#'
#' @examples
#'
estimate_gps <- function(formula,
                         data,
                         method = 'multinom',
                         reference = NULL,
                         by = NULL,
                         missing = NULL,
                         fit.object = FALSE,
                         verbose.output = FALSE,
                         ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  call <- match.call()
  args <- list(...)

  #formula
  if(missing(formula)) {
    chk::abort_chk('The argument `formula` is missing with no default')
  }

  if(!rlang::is_formula(formula, lhs = TRUE)) {
    chk::abort_chk('The argument `formula` has to be a valid R formula with
                   treatment and predictor variables')
  }

  #data



    #if formula_type = subsetting
      # extract data and check if variables in the data
    # if formula type = normal
      # then dataset must be specified
      # check data.frame and variables

  ##

  #method
    # check if method is a one length string
    # check if method in available methods

  #reference
    # check if a one length string
    # has to be in the unique
  #missing

  #by
    #if by is not null, then data has to be provided

  # fit.object + verbose.output
  chk::chk_all(list(fit.object, verbose.output), chk::chk_logical)




  ####################### FITTING ##############################################

  ####################### OUTPUT ###############################################

  #defining the class of the output
}
