#' @title Match the data based on generalized propensity score
#'
#' @description The `match_gps()` function performs sample matching based on
#'   generalized propensity scores (GPS). It utilizes the k-means clustering
#'   algorithm to partition the data into clusters and subsequently matches all
#'   treatment groups within these clusters. This approach ensures efficient and
#'   structured comparisons across treatment levels while accounting for the
#'   propensity score distribution.
#'
#' @param csmatrix An object of class `gps` and/or `csr` representing a data
#'   frame of generalized propensity scores. The first column must be the
#'   treatment variable, with additional attributes describing the calculation
#'   of the common support region and the estimation of generalized propensity
#'   scores. It is crucial that the common support region was calculated using
#'   the `csregion()` function to ensure compatibility.
#' @param caliper A numeric value specifying the caliper width, which defines
#'   the allowable range within which observations can be matched. It is
#'   expressed as a percentage of the standard deviation of the
#'   logit-transformed generalized propensity scores. To perform matching
#'   without a caliper, set this parameter to a very large value. For exact
#'   matching, set `caliper = 0` and enable the `exact` option by setting it to
#'   `TRUE`.
#' @param reference A single string specifying the exact level of the treatment
#'   variable to be used as the reference in the matching process. All other
#'   treatment levels will be matched to this reference level. Ideally, this
#'   should be the control level. If no natural control is present, avoid
#'   selecting a level with extremely low or high covariate or propensity score
#'   values. Instead, choose a level with covariate or propensity score
#'   distributions that are centrally positioned among all treatment groups to
#'   maximize the number of matches.
#' @param combos An optional `data.frame` with two columns. Each row specifies a
#'   single execution of the matching algorithm, where the treatment level in
#'   the second column is matched to the treatment level in the first column
#'   (reference). Matching is performed within clusters calculated by the
#'   k-means algorithm, utilizing generalized propensity scores from all other
#'   treatment levels, excluding the two specified in the `combos` data frame.
#'   The data frame must contain valid strings representing treatment levels,
#'   and its number of rows must be fewer than the total number of unique
#'   pairwise combinations of treatment levels.
#' @param ratio A scalar for the number of matches which should be found. The
#'   default is one-to-one matching.
#' @param replace Logical value indicating whether matching should be done with
#'   replacement. If `FALSE`, the order of matches generally matters. Matches
#'   are found in the same order as the data is sorted. Specifically, the
#'   matches for the first observation will be found first, followed by those
#'   for the second observation, and so on. Matching without replacement is
#'   generally not recommended as it tends to increase bias. However, in cases
#'   where the dataset is large and there are many potential matches, setting
#'   `replace = FALSE` often results in a substantial speedup with negligible or
#'   no bias.
#' @param kmeans_args A list of arguments to pass to [stats::kmeans]. These
#'   arguments must be provided inside a `list()` in the paired `name = value`
#'   format.
#' @param kmeans_cluster An integer specifying the number of clusters to pass to
#'   [stats::kmeans].
#' @param verbose_output a logical flag. If `TRUE` a more verbose version of the
#'   function is run and the output is printed out to the console.
#' @param ... Additional arguments to be passed to the [Matching::Matchby()]
#'   function.
#'
#' @returns A `data.frame` similar to the one provided as the `data` argument in
#'   the [estimate_gps()] function, containing the same columns but only the
#'   observations for which a match was found. The returned object includes two
#'   attributes, accessible with the `attr()` function:
#' * `original_data`: A `data.frame` with the original data returned by the
#'   [csregion()] or [estimate_gps()] function, after the estimation of the csr
#'   and filtering out observations not within the csr.
#'
#' * `matching_filter`: A logical vector indicating which rows from
#'   `original_data` were included in the final matched dataset.
#'
#' @references Michael J. Lopez, Roee Gutman "Estimation of Causal Effects with
#' Multiple Treatments: A Review and New Ideas," Statistical Science, Statist.
#' Sci. 32(3), 432-454, (August 2017)
#'
#' @seealso [estimate_gps()] for the calculation of generalized propensity
#' scores; [Matching::Matchby()] for the documentation of the matching function;
#' [stats::kmeans()] for the documentation of the k-Means algorithm.
#'
#' @examples
#' # Loading the lalonde dataset from the `Matching` package
#' library(Matching)
#' data(lalonde)
#'
#' # Defining the formula used for gps estimation
#' formula_lalonde <- formula(treat ~ age * black * married + educ)
#'
#' # Step 1.) Estimation of the generalized propensity scores
#' gp_scores <- estimate_gps(formula_lalonde,
#'   data = lalonde,
#'   method = "brglm2",
#'   reference = "0",
#'   verbose_output = TRUE
#' )
#'
#' # Step 2.) Defining the common support region
#' gps_csr <- csregion(gp_scores)
#'
#' # Step 3.) Matching the gps
#' matched_lalonde <- match_gps(gps_csr,
#'   caliper = 100, # very big caliper to let
#'   reference = "0", # the gps match freely
#'   replace = FALSE,
#'   kmeans_cluster = 1,
#'   kmeans_args = list(
#'     iter.max = 200,
#'     algorithm = "Forgy"
#'   ),
#'   verbose_output = TRUE
#' )
#'
#' @export
match_gps <- function(csmatrix = NULL,
                      caliper = 0.2,
                      reference = NULL,
                      combos = NULL,
                      ratio = 1,
                      replace = FALSE,
                      kmeans_args = NULL,
                      kmeans_cluster = 5,
                      verbose_output = FALSE,
                      ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  args <- list(...)
  kmeans_args <- kmeans_args %||% list()

  # check and process the csmatrix
  .chk_cond(
    is.null(csmatrix) || missing(csmatrix),
    "The argument `csmatrix` is missing with no default."
  )

  .chk_cond(
    any(c("gps", "data.frame") %nin% class(csmatrix)),
    "The argument `csmatrix` has to be of classes `gps` and
            `data.frame`."
  )

  # Perform the logit transformation and combine with treatment
  logit_matrix <- cbind(treatment = csmatrix[, 1], logit(csmatrix[, -1]))

  # reference
  csmatrix$treatment <- as.factor(csmatrix$treatment)
  ref_list <- .process_ref(csmatrix$treatment,
    ordinal_treat = NULL,
    reference = reference
  )

  args[["treatment"]] <- ref_list[["data.relevel"]]
  reference <- ref_list[["reference"]]

  # kmeans_args, check list, process later
  if (!is.null(kmeans_args)) {
    chk::chk_list(kmeans_args)
  }

  # process combos
  if (is.null(combos)) {
    # generate all possible matches with reference on the left
    combos <- expand.grid(group1 = reference, group2 = colnames(csmatrix)[-1])
    combos <- combos[combos[, 2] != reference, ]
  } else {
    # check if data frame
    .check_df(combos, data_name = "combos")

    # check if two cols
    .chk_cond(
      ncol(combos) != 2,
      "The `combos` data frame must have exactly 2 columns."
    )

    # check if all values in colnames(csmatrix)
    vectorized_combos <- as.character(c(combos[, 1], combos[, 2]))
    .chk_cond(
      any(vectorized_combos %nin% colnames(csmatrix)[-1]),
      "All values in the `combos` table must match the names of the
              unique levels of the treatment variable or the column names of
              the `csmatrix`."
    )

    # check if no repeats (e.g. 1, 1)
    combos_equal <- combos[, 1] == combos[, 2]
    .chk_cond(
      any(combos_equal),
      sprintf(
        "You can not match a group to itself, rows: %s",
        word_list(which(combos_equal))
      )
    )

    # check if all combinations unique
    # Sort each row and convert it into a string for comparison
    row_sorted <- apply(combos, 1, function(row) {
      paste(sort(row),
        collapse = ","
      )
    })

    # Check for duplicates
    duplicates <- duplicated(row_sorted)

    .chk_cond(
      any(duplicates),
      sprintf(
        "You can not check the same combination twice, rows: %s",
        word_list(which(duplicates))
      )
    )
  }

  # define combos
  colnames(combos) <- c("group1", "group2")

  combos[] <- lapply(combos, as.character)
  matches_n <- nrow(combos)

  # ratio
  if (length(ratio) != 1) {
    .chk_cond(
      !.check_vecl(ratio, matches_n, check_numeric = FALSE),
      "The `ratio` argument must be either a single integer or an atomic
              vector with a length equal to the number of rows in the `combos`
              data frame."
    )
  }

  .check_integer(ratio, x_name = "ratio")

  args[["M"]] <- .vectorize(ratio, matches_n)

  # replace
  if (length(replace) != 1) {
    .chk_cond(
      !.check_vecl(replace, matches_n, check_numeric = FALSE),
      "The `replace` argument must be either a single logical flag or an
              atomic vector with a length equal to the number of rows in the
              `combos` data frame."
    )

    .chk_cond(
      !is.logical(replace) || anyNA(replace),
      "All values in the `replace` argument have to be logical flags."
    )
  } else {
    chk::chk_flag(replace)
  }

  args[["replace"]] <- .vectorize(replace, matches_n)

  # process caliper
  .chk_cond(
    is.null(caliper),
    "The `caliper` argument can not be NULL."
  )

  .chk_cond(
    !(.check_vecl(caliper, leng = matches_n, check_numeric = TRUE) ||
      .check_vecl(caliper, leng = 1, check_numeric = TRUE)),
    "The `caliper` argument must be either a single numeric or an atomic
            vector with a length equal to the number of rows in the `combos`
            data frame."
  )

  .chk_cond(
    any(caliper <= 0),
    "The `caliper` argument has to be a positive number."
  )

  caliper <- caliper * stats::sd(as.matrix(logit_matrix[, -1]))

  args[["caliper"]] <- .vectorize(caliper, matches_n)

  # process kmeans_cluster
  .chk_cond(
    is.null(kmeans_cluster),
    "The `kmeans_cluster` argument can not be NULL."
  )

  .chk_cond(
    !(.check_vecl(kmeans_cluster, leng = matches_n, check_numeric = TRUE) ||
      .check_vecl(kmeans_cluster, leng = 1, check_numeric = TRUE)),
    "The `kmeans_cluster` argument must be either a single integer or an
            atomic vector with a length equal to the number of rows in the
            `combos` data frame."
  )

  .check_integer(kmeans_cluster, x_name = "kmeans_cluster")

  .chk_cond(
    any(kmeans_cluster < 1),
    "The `kmeans_cluster` argument must be an integer
            greater than or equal 1."
  )

  kmeans_args[["centers"]] <- .vectorize(kmeans_cluster, matches_n)

  ## deal with algorithm argument
  if (is.null(kmeans_args[["algorithm"]])) {
    kmeans_args["algorithm"] <- "Hartigan-Wong"
  }

  ## processing the kmeans arglist
  kmeans_args <- match_add_args(
    arglist = kmeans_args,
    funlist = stats::kmeans
  )

  ## vectorize the args
  kmeans_args <- lapply(kmeans_args, .vectorize, matches_n)

  ## check if tolerance specified
  tolerance_spec <- is.null(args[["tolerance"]])

  ## processing the matching arglist
  args <- match_add_args(
    arglist = args,
    funlist = Matching::Matchby
  )

  ## vectorize args
  args <- lapply(args, .vectorize, matches_n)
  args <- args[unlist(lapply(args, function(x) !is.null(x)))]

  ## remove tolerance if not specified earlier
  if (tolerance_spec) args <- args[names(args) %nin% "tolerance"]

  ######################## MATCHING ############################################
  ## Steps
  ## 1.) combos --> check and generate
  ## 2.) expand ratio and replace to vectors
  ## 3.) k.cluster argument --> to specify the number of clusters
  ## 4.) actual matching process
  ## 5.) merging the datasets
  ## 6.) output with special class (vector_matched)
  ## 7.) vecotrize kmeans_args and ... args
  ## 8.) refine logic for checking args

  match_results <- list()
  logit_matrix <- cbind(
    ID = seq_len(nrow(logit_matrix)),
    logit_matrix
  )

  for (i in seq_len(matches_n)) {
    ## select only treatment which in current combos loop
    obs_filter <- logit_matrix[, "treatment"] %in% as.vector(combos[i, ])

    # selecting elements from arguments list based on current iteration
    kmeans_args_loop <- lapply(kmeans_args, function(x) x[[i]])
    args_loop <- lapply(args, function(x) x[[i]])

    # selecting columns for kmeans clustering
    cols_kmeans <- colnames(logit_matrix)[-c(1:2)]
    cols_kmeans <- cols_kmeans[cols_kmeans %nin% as.vector(combos[i, ])]

    # selecting from df
    kmeans_args_loop[["x"]] <- logit_matrix[, cols_kmeans]

    # fitting kmeans clusters
    tryCatch(
      {
        verbosely(
          withr::with_preserve_seed(
            k_res <- do.call(
              stats::kmeans,
              kmeans_args_loop
            )
          ),
          verbose = verbose_output
        )
      },
      error = function(e) {
        chk::abort_chk(strwrap(sprintf(
          "There was a problem fitting the kmeans clustering with
        `stats::kmeans()`.\n Error message: (from `stats::kmeans()`) %s",
          conditionMessage(e)
        ), prefix = " ", initial = ""), tidy = FALSE)
      }
    )

    # adding the clusters to matching arguments
    args_loop[["by"]] <- k_res$cluster[obs_filter]

    # selecting columns for matching --> we need only one column of gps!
    cols_matching <- colnames(logit_matrix) == combos[i, 1]

    # selecting data for matching: X argument
    args_loop[["X"]] <- logit_matrix[obs_filter, cols_matching]
    ids_df <- logit_matrix[obs_filter, "ID"]

    # defining Tr args for matching function
    args_loop[["Tr"]] <- ifelse(logit_matrix[obs_filter, 2] == combos[i, 1],
      0,
      1
    )

    # performing the matching
    tryCatch(
      {
        verbosely(
          withr::with_preserve_seed(
            matched <- do.call(
              Matching::Matchby,
              args_loop
            )
          ),
          verbose = verbose_output
        )
      },
      error = function(e) {
        chk::abort_chk(strwrap(sprintf(
          "There was a problem with matching the samples using
        `Matching::Matchby()`.\n Error message:
        (from `Matching::Matchby()`) %s",
          conditionMessage(e)
        ), prefix = " ", initial = ""), tidy = FALSE)
      }
    )

    ## when no matches found, return an error
    .chk_cond(
      all(is.na(matched)),
      "No matches found!"
    )

    res_df <- data.frame(
      control = ids_df[matched$index.control],
      treatment = ids_df[matched$index.treated]
    )

    colnames(res_df) <- paste0("level_", combos[i, 1:2])

    match_results <- append(match_results, list(res_df))
  }

  names(match_results) <- paste0("combos_", seq_len(matches_n))

  ## Works for all combos:
  match_results_outer <- Reduce(
    function(x, y) merge(x, y, all = FALSE),
    match_results
  )

  merged_ids <- unlist(match_results_outer)

  # returning logical vector
  matched_filter <- seq_len(nrow(csmatrix)) %in% merged_ids

  # return original csr data.frame but matched
  csr_data <- attr(csmatrix, "csr_data")
  csr_data <- as.data.frame(csr_data[matched_filter, ])
  attr(csr_data, "matching_filter") <- matched_filter
  if ("csr" %in% class(csmatrix)) {
    attr(csr_data, "original_data") <- attr(csmatrix, "csr_data")
  } else {
    attr(csr_data, "original_data") <- attr(csmatrix, "original_data")
  }

  # Assign class
  class(csr_data) <- c("data.frame", "gps", "csr", "matched")

  return(csr_data)
}
