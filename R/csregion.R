#' @title Filter the data based on common support region
#'
#' @description The `csregion()` function estimates the boundaries of the
#'   rectangular common support region, as defined by Lopez and Gutman (2017),
#'   and filters the matrix of generalized propensity scores based on these
#'   boundaries. The function returns a matrix of observations whose generalized
#'   propensity scores lie within the treatment group-specific boundaries.
#'
#' @param gps_matrix An object of classes `gps` and `data.frame` (e.g., created
#'   by the `estimate_gps()` function). The first column corresponds to the
#'   treatment or grouping variable, while the other columns represent the
#'   treatment assignment probabilities calculated separately for each
#'   hypotetical treatment group. The number of columns should therefore be
#'   equal to the number of unique levels of the treatment variable plus one
#'   (for the treatment variable itself). The number of rows should correspond
#'   to the number of subjects for which generalized propensity scores were
#'   estimated.
#'
#' @return A numeric matrix similar to the one returned by `estimate_gps()`,
#'   but with the number of rows reduced to exclude those observations that do
#'   not fit within the common support region (CSR) boundaries. The returned
#'   object also possesses additional attributes that summarize the calculation
#'   process of the CSR boundaries:
#'  * `filter_matrix` - A logical matrix with the same dimensions as the
#'  gps-part of `gps_matrix`, indicating which treatment assignment
#'  probabilities fall within the CSR boundaries,
#'  * `filter_vector` - A vector indicating whether each observation was kept
#'  (`TRUE`) or removed (`FALSE`), essentially a row-wise
#'  sum of `filter_matrix`,
#'  * `csr_summary` - A summary of the CSR calculation process, including
#'  details of the boundaries and the number of observations filtered.
#'  * `csr_data` - The original dataset used for the estimation of generalized
#'  propensity scores (`original_data` attribute of the `gps` object) filtered
#'  by the `filter_vector`
#'
#' @examples
#' # We could estimate simples generalized propensity scores for the `iris`
#' # dataset
#' gps <- estimate_gps(Species ~ Sepal.Length, data = iris)
#'
#' # And then define the common support region boundaries using `csregion()`
#' gps_csr <- csregion(gps)
#'
#' # The additional information of the CSR-calculation process are
#' # accessible through the attributes described in the `*Value*` section
#' attr(gps_csr, "filter_matrix")
#' attr(gps_csr, "csr_summary")
#' attr(gps_csr, "csr_data")
#'
#' @export
csregion <- function(gps_matrix) {
  csr_data <- attr(gps_matrix, "original_data")

  .chk_cond(
    "gps" %nin% class(gps_matrix),
    "The `gps_matrix` argument must be of class `gps`."
  )

  ## Calculating the csr_low
  csr_low <- apply(
    stats::aggregate(. ~ treatment,
      data = gps_matrix,
      FUN = function(x) min(x)
    )[, -1],
    2,
    max
  )

  ## Calculating the csr_high
  csr_high <- apply(
    stats::aggregate(. ~ treatment,
      data = gps_matrix,
      FUN = function(x) max(x)
    )[, -1],
    2,
    min
  )

  ## filter out the unvalid observations
  filter_matrix <- mapply(function(df_col, vec_low, vec_high) {
    vec_low < df_col & vec_high > df_col
  }, gps_matrix[, 2:ncol(gps_matrix)], csr_low, csr_high)

  ## summarizing the logical matrix to a subset vector
  filter_vector <- apply(filter_matrix, 1, all)

  ## defining the number of negatives
  n_negative <- sum(!filter_vector)
  n_negative_matrix <- colSums(!filter_matrix)

  ## subsetting the gps_matrix
  gps_matrix <- subset(gps_matrix, filter_vector)

  ## drop unused levels of treatment variable
  gps_matrix[, "treatment"] <- droplevels(gps_matrix[, "treatment"])

  ## detect low number of observations in groups and print a warning
  if (any(table(gps_matrix[, "treatment"]) < 20)) {
    chk::wrn(strwrap("Some groups have fewer than 20 observations, which may
    impact the performance of the matching process. Consider using
    `replace = TRUE`in `match_gps()` to address this.",
      prefix = " ", initial = ""
    ))
  }

  ## assembling the print results list
  csr_summary <- data.frame(
    treatment = colnames(gps_matrix)[2:ncol(gps_matrix)],
    csr_low = csr_low,
    csr_high = csr_high,
    n_negative_matrix = n_negative_matrix
  )

  res <- list(
    data_csr = csr_summary,
    excluded = n_negative
  )

  ## Setting a new class for the results
  csres <- structure(res, class = "csres")

  ## print custom output
  show_csres(csres)

  ## adding attributes to csr output
  attr(gps_matrix, "filter_matrix") <- filter_matrix
  attr(gps_matrix, "filter_vector") <- filter_vector
  attr(gps_matrix, "csr_summary") <- csr_summary
  attr(gps_matrix, "csr_data") <- csr_data[filter_vector, ]

  # Assign attributes and class
  class(gps_matrix) <- c("data.frame", "gps", "csr")

  # return the gps_matrix
  return(invisible(gps_matrix))
}

# Function to show the contents of an object of class 'csres'
show_csres <- function(object) {
  # Function to print a data frame in a table-like format
  print_table <- function(df, colnames_df) {
    # Print column headers
    cat(paste(sprintf("%-15s", colnames_df), collapse = " | "), "\n")
    cat(paste(rep("-", 17 * ncol(df)), collapse = ""), "\n")

    # Print each row of the data frame
    apply(df, 1, function(row) {
      cat(paste(sprintf("%-15s", row), collapse = " | "), "\n")
    })
  }

  # Print header and content
  cat("\n")
  cat("Rectangular CSR Borders Evaluation", "\n")
  cat("==================================\n\n")

  print_table(object$data_csr,
    colnames = c(
      "Treatment", "Lower CSR limit",
      "Upper CSR limit", "Number excluded"
    )
  )

  cat("\n")
  cat("===================================================\n")
  cat("The total number of excluded observations is:\t", object$excluded, "\n")
  cat(strwrap("Note: You can view the summary of the
              CSR calculation using the `attr()` function.",
    prefix = " ", initial = ""
  ))
}
