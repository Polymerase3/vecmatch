#--check names of variables-----------------------------------------------------
## Returns TRUE if object is string character or name of given length
.val_name_q <- function(x, len) {
  return((is.symbol(x) || is.character(x)) && length(eval(x)) == len)
}

#--check if name in the data----------------------------------------------------
.check_name <- function(data, namlist) {
  convert <- unlist(lapply(namlist, function(x) !is.character(x)))
  if(length(convert) != 0) namlist[convert] <- lapply(namlist[convert], as.character)

  which <- unlist(namlist) %in% names(data)
  if(!all(which)) return(unlist(namlist)[!which])
}

#--check if the name is a symbol or a string------------------------------------
.chk_symbol <- function(...) {
  tryCatch(
    {
      symbols <- rlang::ensyms(..., .ignore_null = 'all')
      return(symbols)
    }, error = function() {
      chk::abort_chk('The provided column names are not symbols or text strings!')
    })
}


