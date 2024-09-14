#' Title
#'
#' @return
#' @export
#'
#' @examples
raincloud <- function(data = NULL,
                      y = NULL,
                      group = NULL,
                      facet = NULL,
                      significance = FALSE,
                      limits = NULL,
                      jitter = NULL,
                      alpha = 0.4,
                      save = FALSE) {

  ############################INPUT CHECKING####################################
  #--check data frame-----------------------------------------------------------
  ## Must be an object of class data frame with at least one numeric column
  if(is.null(data) || !inherits(data, 'data.frame')) {
    chk::abort_chk('Argument `data` must be an object of class `data.frame`')
  }

  if(length(data) == 0) {
    chk::abort_chk('The provided data frame is empty')
  }

  if(length(data) == 1 && !is.numeric(data[, 1])) {
    chk::abort_chk('The provided data is not numeric')
  }

  #--check y, group and facet---------------------------------------------------

  ## Check if the provided names are valid names + convert to
  symlist <- list(y = substitute(y),
                  group = substitute(group),
                  facet = substitute(facet))
  symlist <- .conv_nam(symlist)

  ## Check if y exists
  if(is.null(symlist[[1]])) chk::abort_chk('The argument `y` is missing
                                           with no default!')

  ## Check if there are in the dataframe
  nonames <- .check_name(data, symlist)

  if(length(nonames) != 0) {
    chk::abort_chk(paste0('The following colnames are not in the
                          provided data frame: ',
                          paste(nonames, collapse = ', ')
                          ))
  }

  ## Check if significance is logical
  chk::chk_logical(significance)

  ## Check if limits is a numeric vector of length 2
  if(!is.null(limits)) {
    if(!.check_vecl(limits, leng = 2)) {
      chk::abort_chk('The `limits` argument should be a numeric vector of length
                     2: c(`min`, `max`)')
    }
  }

  ## Check range for jitter
  if(!is.null(jitter)) chk::chk_range(jitter, range = c(0, 1))

  ## Check range for alpha
  chk::chk_range(alpha, range = c(0, 1))

  ## Check logical for save
  chk::chk_logical(save)

  ####################### DATA PROCESSING ######################################
  # assure y is numeric and convert facet, group to factors
  mapply(.conv_data,
         type = list('numeric', 'factor', 'factor'),
         varname = symlist,
         MoreArgs = list(data = data,
                         env = environment()))

  ####################### PLOTTING #############################################
}
