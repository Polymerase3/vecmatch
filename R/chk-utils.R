#--check if name in the data----------------------------------------------------
.check_name <- function(data, namlist) {
  convert <- unlist(lapply(namlist, function(x) !is.character(x)))
  if (length(convert) != 0) {
    namlist[convert] <- lapply(
      namlist[convert],
      as.character
    )
  }

  which <- unlist(namlist) %in% names(data)
  if (!all(which)) {
    return(unlist(namlist)[!which])
  }
}

#--check if the object is a numeric vector of given length----------------------
.check_vecl <- function(x, leng) {
  ret <- all(
    is.atomic(x) && !is.matrix(x) && !is.array(x), # checks vector
    length(x) == leng, # checks length
    is.numeric(x) # checks numeric
  )
  return(ret)
}

#--convert list with names to character and leave NULLs-------------------------
.conv_nam <- function(namlist, check_valid = TRUE) {
  # if NULL then leave
  # if not null then quote
  .conv <- function(x) {
    if (!is.null(x)) {
      x <- rlang::quo_name(x)
      if (check_valid == TRUE) chk::chk_valid_name(x)
      return(x)
    }
  }
  res <- lapply(namlist, .conv)
  return(res)
}

#--convert data with error message----------------------------------------------
.conv_data <- function(data, type, varname, env) {
  if (is.null(varname)) invisible(return())

  column <- paste0('data[, "', varname, '"]')
  cond <- paste0("!is.", type, "(", column, ")")
  conv_call <- paste0(column, "<- as.", type, "(", column, ")")
  abort_call <- quote(chk::abort_chk(paste0(
    "The variable `", varname,
    "` cannot be converted to the type ",
    type, "."
  )))

  if (type == "factor") {
    def_type <- (is.numeric(eval(parse(text = column))) ||
      is.integer(eval(parse(text = column))) ||
      is.factor(eval(parse(text = column))) ||
      is.character(eval(parse(text = column))) ||
      is.logical(eval(parse(text = column))))
    if (!def_type) eval(abort_call)

    length_unique <- length(unique(eval(parse(text = column))))

    if (length_unique > 10 && def_type) {
      chk::wrn(paste0(
        "The variable ", varname,
        "has more tha 10 unique values. ",
        "Are you sure you want to proceed?"
      ))
    }
  }

  env$conv_call <- conv_call
  if (eval(parse(text = cond))) {
    tryCatch(
      {
        invisible(with(env, {
          eval(parse(text = conv_call))
        }))
      },
      warning = function(w) {
        eval(abort_call)
      }
    )
  }
}

#--check if filename has given extension
.check_extension <- function(name, x_name, ext_vec) {
  ext <- substr(name, nchar(name) - 3, nchar(name))
  cond <- ext %in% ext_vec
  if (!cond) {
    chk::abort_chk(sprintf('The argument `%s` is allowed to have the following extensions: %s',
                           x_name, word_list(add_quotes(ext_vec))))
  }
}

#--check data.frame-------------------------------------------------------------
.check_df <- function(data) {
  if (is.null(data) || !inherits(data, "data.frame")) {
    chk::abort_chk("Argument `data` must be an object of class `data.frame`")
  }

  if (length(data) == 0) {
    chk::abort_chk("The provided data frame is empty")
  }
}

#--check gps methods------------------------------------------------------------
.check_method <- function(string) {
  # if (missing(string) || is.null(string)) invisible(return(NULL))

  if (!(is.character(string) && length(string) == 1L && !anyNA(string))) {
    chk::abort_chk("The argument `method` must be a single string of length 1")
  }

  if (!(string %in% names(.gps_methods))) {
    chk::abort_chk(sprintf(
      "The argument `method` has to be one from: %s",
      word_list(add_quotes(names(.gps_methods)))
    ))
  }
}

#--process link function--------------------------------------------------------
