#--check custom condition-------------------------------------------------------
.chk_cond <- function(condition, error_message, error = TRUE, ...) {
  if (condition && error) {
    chk::abort_chk(
      strwrap(error_message, prefix = " ", initial = ""),
      ...
    )
  } else if (condition && !error) {
    chk::wrn(
      strwrap(error_message, prefix = " ", initial = ""),
      ...
    )
  }
}

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
.check_vecl <- function(x, leng = NULL, check_numeric = TRUE) {
  ret <- all(
    is.atomic(x) && !is.matrix(x) && !is.array(x), # checks vector
    if (!is.null(leng)) length(x) == leng, # checks length
    if (check_numeric) is.numeric(x) # checks numeric
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
  abort_call <- quote(
    chk::abort_chk(
      strwrap(
        sprintf(
          "The variable `%s` cannot be converted to the type %s",
          varname, type
        ),
        prefix = " ", initial = ""
      )
    )
  )

  if (type == "factor") {
    def_type <- (is.numeric(eval(parse(text = column))) ||
      is.integer(eval(parse(text = column))) ||
      is.factor(eval(parse(text = column))) ||
      is.character(eval(parse(text = column))) ||
      is.logical(eval(parse(text = column))))
    if (!def_type) eval(abort_call)

    length_unique <- length(unique(eval(parse(text = column))))

    if (length_unique > 10 && def_type) {
      chk::wrn(strwrap(sprintf("The variable %s has more tha 10 unique values.
                               Are you sure you want to proceed?", varname),
        prefix = " ", initial = ""
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
    chk::abort_chk(strwrap(sprintf(
      "The argument `%s` is allowed to have the following extensions: %s",
      x_name, word_list(add_quotes(ext_vec))
    ), prefix = " ", initial = ""))
  }
}

#--check data.frame-------------------------------------------------------------
.check_df <- function(data, data_name = "data") {
  .chk_cond(
    is.null(data) || !inherits(data, "data.frame") ||
      length(class(data)) > 1,
    sprintf("Argument `%s` must be an object of
                    single class `data.frame`", data_name)
  )

  .chk_cond(length(data) == 0, "The provided data frame is empty")
}

#--check gps methods------------------------------------------------------------
.check_method <- function(string) {
  .chk_cond(
    !(is.character(string) && length(string) == 1L && !anyNA(string)),
    "The argument `method` must be a single string of length 1"
  )

  .chk_cond(
    !(string %in% names(.gps_methods)),
    sprintf(
      "The argument `method` has to be one from: %s",
      word_list(add_quotes(names(.gps_methods)))
    )
  )
}

## --check gps matrix-----------------------------------------------------------
.check_gps_matrix <- function(gps_matrix) {
  # check if treatment variable present and if its first column of all
  .chk_cond(
    colnames(gps_matrix)[1] != "treatment",
    "The first column in an object of class `gps` must be named
            `treatment`."
  )

  # check if all levels of treatment are in colnames
  .chk_cond(
    nunique(gps_matrix[, 1]) != ncol(gps_matrix) - 1,
    "The `gps` object must have a number of columns equal to the unique
            levels of the treatment variable plus one (to include the treatment
            variable itself)."
  )

  .chk_cond(
    any(unique(gps_matrix$treatment) %nin% colnames(gps_matrix)),
    "The columns of the `gps` object must be named to match the unique
            levels of the treatment variable."
  )

  .chk_cond(
    any(is.na(gps_matrix)),
    "The object of class `gps` can not contain any NA's"
  )

  if (any(round(rowSums(gps_matrix[, -1]), 5) != 1)) {
    which_row <- which(round(rowSums(gps_matrix[, -1]), 5) != 1)

    chk::abort_chk(strwrap(sprintf(
      "The row-wise sum of probabilities across all columns must equal 1.
      IDs of rows where this condition is not met: %s",
      word_list(which_row)
    ), prefix = " ", initial = ""))
  }
}

## --check integer--------------------------------------------------------------
.check_integer <- function(x, x_name = NULL) {
  coerce_integer <- suppressWarnings(as.integer(x))

  if (!any(is.na(coerce_integer))) {
    is_integer <- all.equal(x, coerce_integer, giveErr = TRUE)

    if (!is.null(attr(is_integer, "err", exact = TRUE))) is_integer <- FALSE
  } else {
    is_integer <- FALSE
  }

  .chk_cond(!is_integer, sprintf(
    "The argument `%s` has to be an integer",
    x_name
  ))
}

## --match discrete args--------------------------------------------------------
.match_discrete_args <- function(x, choices, x_name) {
  tryCatch(
    {
      unlist(lapply(x, match.arg, choices = choices))
    },
    error = function(e) {
      chk::abort_chk(strwrap(sprintf(
        "The `%s` argument can only have following
                                   values: %s", x_name,
        word_list(
          add_quotes(
            choices
          )
        )
      )))
    }
  )
}

## --check null and replace with default----------------------------------------
.chk_null_default <- function(x, x_name, method, default) {
  .chk_cond(is.null(x),
    error = FALSE,
    sprintf(
      "The `%s` argument for the method %s was not provided
                    and will default to `%s`.",
      x_name,
      add_quotes(method),
      default
    )
  )

  if (is.null(x)) x <- default

  return(x)
}

.chk_vararg_length <- function(x, x_name,
                               check_numeric = FALSE,
                               type_n = "integer",
                               matches_n) {
  .chk_cond(
    !(.check_vecl(x, matches_n, check_numeric = check_numeric) ||
      .check_vecl(x, leng = 1, check_numeric = check_numeric)),
    sprintf("The `%s` argument must be either a single %s or an atomic
            vector with a length equal to the number of rows in the `combos`
            data frame.", x_name, type_n)
  )
}
