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
                      limits = NULL,
                      jitter = NULL,
                      alpha = NULL,
                      legend.position = NULL,
                      kernel.method = NULL,
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
  tryCatch(
    {
      symlist <- rlang::ensyms(y, facet, group, .ignore_null = 'all')
    }, error = function(e) {
      chk::abort_chk('The provided column names are not symbols or text strings!')
    })

  ## Check if there are in the dataframe
  nonames <- .check_name(data, symlist)

  if(length(nonames) != 0) {
    chk::abort_chk(paste0('The following colnames are not in the
                          provided data frame: ',
                          paste(nonames, collapse = ', ')
                          ))
  }

}
