#--check if name in the data----------------------------------------------------
.check_name <- function(data, namlist) {
  convert <- unlist(lapply(namlist, function(x) !is.character(x)))
  if(length(convert) != 0) namlist[convert] <- lapply(namlist[convert],
                                                      as.character)

  which <- unlist(namlist) %in% names(data)
  if(!all(which)) return(unlist(namlist)[!which])
}

#--check if the object is a numeric vector of given length----------------------
.check_vecl <- function(x, leng) {
  ret <- all(
    is.atomic(x) && !is.matrix(x) && !is.array(x), # checks vector
    length(x) == leng,                             # checks length
    is.numeric(x)                                  # checks numeric
  )
  return(ret)
}

#--convert list with names to character and leave NULLs-------------------------
.conv_nam <- function(namlist, check_valid = TRUE) {
  # if NULL then leave
  # if not null then quote
  .conv <- function(x) {
    if(!is.null(x)) {
      x <- rlang::quo_name(x)
      if(check_valid == TRUE) chk::chk_valid_name(x)
      return(x)
    }
  }
  res <- lapply(namlist, .conv)
  return(res)
}

#--convert data with error message----------------------------------------------

