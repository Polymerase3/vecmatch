#--check if name in the data----------------------------------------------------
.check_name <- function(data, namlist) {
  convert <- unlist(lapply(namlist, function(x) !is.character(x)))
  if(length(convert) != 0) namlist[convert] <- lapply(namlist[convert], as.character)

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

