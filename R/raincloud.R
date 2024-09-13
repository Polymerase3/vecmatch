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


  ## Check if the provided names are symbols or strings
  y <- rlang::quo_name(rlang::enquo(y))
  group <- rlang::quo_name(rlang::enquo(group))
  facet <- rlang::quo_name(rlang::enquo(facet))

  print(y)
  chk::chk_valid_name()

  ## make named list from names
  ## error if y == 'NULL'
  ## write a function to handle 'NULLs'


  tryCatch(
    {
      symlist <- rlang::ensyms(y, facet, group, .ignore_null = 'all')
    }, error = function(e) {
      chk::abort_chk('The provided column names are not symbols or text strings!')
    })

  ## Check if y exists
  if(is.null(rlang::ensym(y))) chk::abort_chk('The argument `y` is missing with no default!')



  #names(symlist)[1] <- 'y'

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

  ####################### DATA PROCESSING#######################################

}
