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
#' @param method A single string specifying the matching method to use. The
#'   default is `"nnm"`, which applies the k-nearest neighbors matching
#'   algorithm. See the Details section for a full list of available methods.
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
#' @param ratio A scalar for the number of matches which should be found for
#'   each control observation. The default is one-to-one matching. Only
#'   available for the methods `"nnm"` and `"pairopt"`.
#' @param replace Logical value indicating whether matching should be done with
#'   replacement. If `FALSE`, the order of matches generally matters. Matches
#'   are found in the same order as the data is sorted. Specifically, the
#'   matches for the first observation will be found first, followed by those
#'   for the second observation, and so on. Matching without replacement is
#'   generally not recommended as it tends to increase bias. However, in cases
#'   where the dataset is large and there are many potential matches, setting
#'   `replace = FALSE` often results in a substantial speedup with negligible or
#'   no bias. Only available for the method `"nnm"`
#' @param order A string specifying the order in which logit-transformed GPS
#'   values are sorted before matching. The available options are:
#'  * `"desc"` – sorts GPS values from highest to lowest (default).
#'  * `"asc"` – sorts GPS values from lowest to highest.
#'  * `"original"` – preserves the original order of GPS values.
#'  * `"random"` – randomly shuffles GPS values. To generate different random
#'  orders, set a seed using [set.seed()].
#' @param ties A logical flag indicating how tied matches should be handled.
#'   Available only for the `"nnm"` method, with a default value of `FALSE` (all
#'   tied matches are included in the final dataset, but only unique
#'   observations are retained). For more details, see the `ties` argument in
#'   [Matching::Matchby()].
#' @param min_controls The minimum number of treatment observations that should
#'   be matched to each control observation. Available only for the `"fullopt"`
#'   method. For more details, see the `min.controls` argument in
#'   [optmatch::fullmatch()].
#' @param max_controls The maximum number of treatment observations that can be
#'   matched to each control observation. Available only for the `"fullopt"`
#'   method. For more details, see the `max.controls` argument in
#'   [optmatch::fullmatch()].
#' @param kmeans_args A list of arguments to pass to [stats::kmeans]. These
#'   arguments must be provided inside a `list()` in the paired `name = value`
#'   format.
#' @param kmeans_cluster An integer specifying the number of clusters to pass to
#'   [stats::kmeans].
#' @param verbose_output a logical flag. If `TRUE` a more verbose version of the
#'   function is run and the output is printed out to the console.
#' @param ... Additional arguments to be passed to the matching
#'   function.
#'
#' @details Propensity score matching can be performed using various matching
#'   algorithms. Lopez and Gutman (2017) do not explicitly specify the matching
#'   algorithm used, but it is assumed they applied the commonly used k-nearest
#'   neighbors matching algorithm, implemented as `method = "nnm"`. However,
#'   this algorithm can sometimes be challenging to use, especially when
#'   treatment and control groups have unequal sizes. When `replace = FALSE`,
#'   the number of matches is strictly limited by the smaller group, and even
#'   with `replace = TRUE`, the results may not always be satisfactory. To
#'   address these limitations, we have implemented an additional matching
#'   algorithm to maximize the number of matched observations within a dataset.
#'
#'   The available matching methods are:
#'
#'   * `"nnm"` – classic k-nearest neighbors matching, implemented using
#'   [Matching::Matchby()]. The tunable parameters in `match_gps()` are
#'   `caliper`, `ratio`, `replace`, `order`, and `ties`. Additional arguments
#'   can be passed to [Matching::Matchby()] via the `...` argument.
#'   * `"fullopt"` – optimal full matching algorithm, implemented with
#'   [optmatch::fullmatch()]. This method calculates a discrepancy matrix to
#'   identify all possible matches, often optimizing the percentage of matched
#'   observations. The available tuning parameters are `caliper`,
#'   `min_controls`, and `max_controls`.
#'   * `"pairmatch"` – optimal 1:1 and 1:k matching algorithm, implemented using
#'   [optmatch::pairmatch()], which is actually a wrapper around
#'   [optmatch::fullmatch()]. Like `"fullopt"`, this method calculates a
#'   discrepancy matrix and finds matches that minimize its sum. The available
#'   tuning parameters are `caliper` and `ratio`.

#'
#' @return A `data.frame` similar to the one provided as the `data` argument in
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
#'   scores; [MatchIt::matchit()], [optmatch::fullmatch()] and
#'   [optmatch::pairmatch()] for the documentation of the matching functions;
#'   [stats::kmeans()] for the documentation of the k-Means algorithm.
#'
#' @examples
#' # Defining the formula used for gps estimation
#' formula_cancer <- formula(status ~ age + sex)
#'
#' # Step 1.) Estimation of the generalized propensity scores
#' gp_scores <- estimate_gps(formula_cancer,
#'   data = cancer,
#'   method = "multinom",
#'   reference = "control",
#'   verbose_output = TRUE
#' )
#'
#' # Step 2.) Defining the common support region
#' gps_csr <- csregion(gp_scores)
#'
#' # Step 3.) Matching the gps
#' matched_cancer <- match_gps(gps_csr,
#'   caliper = 0.25,
#'   reference = "control",
#'   method = "fullopt",
#'   kmeans_cluster = 2,
#'   kmeans_args = list(
#'     iter.max = 200,
#'     algorithm = "Forgy"
#'   ),
#'   verbose_output = TRUE
#' )
#'
#' @export
match_gps <- function(csmatrix = NULL,
                      method = "nnm",
                      caliper = 0.2,
                      reference = NULL,
                      ratio = NULL,
                      replace = NULL,
                      order = NULL,
                      ties = NULL,
                      min_controls = NULL,
                      max_controls = NULL,
                      kmeans_args = NULL,
                      kmeans_cluster = 5,
                      verbose_output = FALSE,
                      ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################

  # capture arguments
  all_args <- as.list(environment())
  ellipsis_args <- list(...)
  kmeans_args <- kmeans_args %||% list()
  args <- list()

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

  # defining the main data.frame
  data_name <- .match_methods[[method]]$data_name
  args[[data_name]] <- logit_matrix

  # Reset the rows counter on the logit_matrix
  rownames(args[[data_name]]) <- NULL

  # Control the method argument
  .chk_cond(
    method %nin% names(.match_methods),
    sprintf(
      "The argument `method` is only allowed to have following
             values: %s",
      word_list(names(.match_methods), quotes = TRUE)
    )
  )

  # Process the variable args for different methods
  varargs <- c(
    "caliper", "ratio", "replace", "order",
    "ties", "min_controls", "max_controls"
  )

  # define which varags are possible for given methods (both named args and ...)
  allowed_varargs <- unlist(
    lapply(
      .match_methods[[method]]$args_check_fun,
      function(x) {
        forms <- formals(x)
        forms[unlist(lapply(forms, is.null))] <- NA
        names(forms)
      }
    )
  )

  allowed_varargs <- c(
    .match_methods[[method]]$allowed_args,
    unique(allowed_varargs)
  )

  # check which varargs were specified (both names and ...)
  specified_varargs <- Filter(Negate(is.null), all_args)
  varargs_filter <- names(specified_varargs) %in% varargs
  specified_varargs <- c(
    specified_varargs[varargs_filter],
    ellipsis_args
  )

  # filter out only allowed args from the specified ones
  checked_varargs_filter <- names(specified_varargs) %in% allowed_varargs

  # if not allowed varargs were specified:
  .chk_cond(
    !all(checked_varargs_filter),
    error = FALSE,
    sprintf(
      "Following specified arguments are not
                    allowed for the method %s and will be ignored: %s",
      add_quotes(method),
      word_list(
        names(specified_varargs[!checked_varargs_filter])
      )
    )
  )

  # reference processing
  csmatrix$treatment <- as.factor(csmatrix$treatment)
  ref_list <- .process_ref(csmatrix$treatment,
    ordinal_treat = NULL,
    reference = reference
  )

  if (.match_methods[[method]]$treat_var) {
    args[["Tr"]] <- ref_list[["data.relevel"]]
  }

  reference <- ref_list[["reference"]]

  # kmeans_args, check list, process later
  if (!is.null(kmeans_args)) {
    chk::chk_list(kmeans_args)
  }

  # process combos
  if ("combos" %nin% names(ellipsis_args)) {
    # generate all possible matches with reference on the left
    combos <- expand.grid(group1 = reference, group2 = colnames(csmatrix)[-1])
    combos <- combos[combos[, 2] != reference, ]
  } else {
    # extracting combos
    combos <- ellipsis_args[["combos"]]

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

  ## varargs processing section ================================================
  # ratio
  if ("ratio" %in% allowed_varargs) {
    ratio <- .chk_null_default(ratio, "ratio", method, 1)

    .chk_vararg_length(ratio, "ratio",
      type_n = "integer",
      matches_n = matches_n
    )
    .check_integer(ratio, x_name = "ratio")

    args[["M"]] <- .vectorize(ratio, matches_n)
  }

  # replace
  if ("replace" %in% allowed_varargs) {
    replace <- .chk_null_default(replace, "replace", method, FALSE)

    .chk_vararg_length(replace, "replace",
      type_n = "logical flag",
      matches_n = matches_n
    )

    .chk_cond(
      !is.logical(replace) || anyNA(replace),
      "All values in the `replace` argument have to be logical flags."
    )

    args[["replace"]] <- .vectorize(replace, matches_n)
  }

  # process caliper
  if ("caliper" %in% allowed_varargs) {
    caliper <- .chk_null_default(caliper, "caliper", method, 0.25)

    .chk_vararg_length(caliper, "caliper",
      TRUE,
      type_n = "numeric",
      matches_n
    )

    .chk_cond(
      any(caliper <= 0),
      "The `caliper` argument has to be a positive number."
    )

    caliper <- caliper * stats::sd(as.matrix(args[[data_name]][, reference]))

    args[["caliper"]] <- .vectorize(caliper, matches_n)
  }

  # define the matching formula based on combos
  if (.match_methods[[method]]$formula_necessary) {
    matching_formulas <- apply(
      combos,
      1,
      function(row) {
        paste0("treatment", " ~ ", row[1])
      }
    )
    args[["x"]] <- unlist(matching_formulas)
  }

  # processing the ties argument
  if ("ties" %in% allowed_varargs) {
    ties <- .chk_null_default(ties, "ties", method, TRUE)

    .chk_vararg_length(ties, "ties",
      type_n = "logical flag",
      matches_n = matches_n
    )

    .chk_cond(
      !is.logical(ties) || anyNA(ties),
      "All values in the `ties` argument have to be logical flags."
    )

    args[["ties"]] <- .vectorize(ties, matches_n)
  }

  # processing the order argument
  if ("order" %in% allowed_varargs) {
    order <- .chk_null_default(order, "order", method, "desc")

    .chk_vararg_length(order, "order",
      type_n = "text string",
      matches_n = matches_n
    )

    allowed_orders <- c("desc", "asc", "original", "random")

    .chk_cond(
      order %nin% allowed_orders,
      sprintf(
        "The `order` argument has to be one of the following values: %s",
        word_list(allowed_orders, quotes = TRUE)
      )
    )

    sort_before_matching <- .vectorize(order, matches_n)
  }

  # processing the min_controls
  if ("min_controls" %in% allowed_varargs) {
    min_controls <- .chk_null_default(min_controls, "min_controls", method, 0)

    .chk_vararg_length(min_controls, "min_controls",
      TRUE,
      type_n = "numeric",
      matches_n
    )

    .chk_cond(
      any(min_controls < 0),
      "The `min_controls` argument has to be a non-negative finite number."
    )

    args[["min.controls"]] <- .vectorize(min_controls, matches_n)
  }

  # processing the max_controls
  if ("max_controls" %in% allowed_varargs) {
    max_controls <- .chk_null_default(max_controls, "max_controls", method, Inf)

    .chk_vararg_length(max_controls, "max_controls",
      TRUE,
      type_n = "numeric",
      matches_n
    )

    .chk_cond(
      any(max_controls <= 0),
      "The `max_controls` argument has to be a non-negative number."
    )

    args[["max_controls"]] <- .vectorize(caliper, matches_n)
  }

  # process kmeans_cluster
  .chk_cond(
    is.null(kmeans_cluster),
    "The `kmeans_cluster` argument can not be NULL."
  )

  .chk_vararg_length(kmeans_cluster, "kmeans_cluster",
    TRUE,
    type_n = "integer",
    matches_n
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

  ## change ratio to controls in "pairopt"
  if (method == "pairopt") names(args)[names(args) == "M"] <- "controls"

  ## processing the matching arglist
  args <- match_add_args(
    arglist = args,
    funlist = .match_methods[[method]]$args_check_fun
  )

  if (method == "nnm") {
    args <- args[names(args) %nin% c("Z", "V", "tolerance")]
  }

  args <- Filter(function(x) !all(is.na(x)), args)

  ## vectorize args
  args <- lapply(args, .vectorize, matches_n)

  ######################## MATCHING ############################################
  match_results <- list()

  for (i in seq_len(matches_n)) {
    ## select only treatment which in current combos loop
    obs_filter <- args[[data_name]][, "treatment"] %in% as.vector(combos[i, ])

    # selecting elements from arguments list based on current iteration
    kmeans_args_loop <- lapply(kmeans_args, function(x) x[[i]])

    # doing the same for args, but without data args
    args_length <- lengths(args)
    args_length[names(args_length) %in% c("Tr", "data", "X")] <- 0

    args_loop <- mapply(
      function(x, y) {
        if (x == matches_n) {
          y[[i]]
        } else {
          y
        }
      }, args_length,
      args,
      SIMPLIFY = FALSE
    )

    # selecting columns for kmeans clustering
    cols_kmeans <- colnames(args[[data_name]])[-1]
    cols_kmeans <- cols_kmeans[cols_kmeans %nin% as.vector(combos[i, ])]

    # selecting from df
    kmeans_args_loop[["x"]] <- args[[data_name]][, cols_kmeans]

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

    ############################ FIXING AREA ###################################
    # adding the clusters to matching arguments
    kmeans_strata <- as.factor(k_res$cluster[obs_filter])

    # selecting columns for matching --> we need only one column of gps!
    cols_matching <- colnames(args[[data_name]]) %in% c(
      "treatment",
      combos[i, 1]
    )

    # arg processing for Matching::Matchby

    if (method == "nnm") {
      # filtering the data frame for matching
      args_loop[[data_name]] <- args_loop[[data_name]][, -1, drop = FALSE]
      args_loop[[data_name]] <- args_loop[[data_name]][
        obs_filter,
        cols_matching
      ]

      # defining data order
      order_data <- .ordering_func(args_loop[[data_name]][, reference],
        order = sort_before_matching[i]
      )
      order_original <- seq_len(nrow(args_loop[[data_name]]))
      order_original <- order_original[order_data]

      # reordering the data
      args_loop[[data_name]] <- args_loop[[data_name]][order_data, ]

      # adding the clusters to matching arguments
      args_loop[["by"]] <- k_res$cluster[obs_filter]
      args_loop[["by"]] <- args_loop[["by"]][order_data]

      # defining Tr args for matching function
      args_loop[["Tr"]] <- args_loop[["Tr"]][obs_filter]
      args_loop[["Tr"]] <- ifelse(args_loop[["Tr"]] == combos[i, 1],
        1,
        0
      )
      args_loop[["Tr"]] <- args_loop[["Tr"]][order_data]
    } else if (method %in% c("fullopt", "pairopt")) {
      # selecting observations and adding "by" argument from kmeans
      args_loop[[data_name]] <- args[[data_name]][obs_filter, cols_matching]
      args_loop[[data_name]] <- cbind(args_loop[[data_name]],
        kmeans_strata = kmeans_strata
      )

      args_loop[[data_name]] <- as.data.frame(args_loop[[data_name]])

      # defining the order
      order_data <- .ordering_func(args_loop[[data_name]][, reference],
        order = order
      )
      order_original <- seq_len(nrow(args_loop[[data_name]]))
      order_original <- order_original[order_data]

      # reordering the data
      args_loop[[data_name]] <- args_loop[[data_name]][order_data, ]

      # converting treatment to binary
      args_loop[[data_name]][, "treatment"] <- ifelse(
        args_loop[[data_name]][, "treatment"] == combos[i, 1],
        0,
        1
      )

      # converting the formula to add kmeans stratas
      args_loop[["x"]] <- paste0(
        args_loop[["x"]],
        " + optmatch::strata(kmeans_strata)"
      )
      args_loop[["x"]] <- stats::as.formula(args_loop[["x"]])

      # performing the prematching =============================================
      tryCatch(
        {
          verbosely(
            withr::with_preserve_seed(
              prematched <- do.call(
                .match_methods[[method]]$args_check_fun[[1]],
                args_loop
              )
            ),
            verbose = verbose_output
          )
        },
        error = function(e) {
          chk::abort_chk(strwrap(sprintf(
            "There was a problem with matching the samples using
        `%s`.\n Error message:
        (from `%s`) %s",
            deparse(substitute(.match_methods[[method]]$args_check_fun[[1]])),
            deparse(substitute(.match_methods[[method]]$args_check_fun[[1]])),
            conditionMessage(e)
          ), prefix = " ", initial = ""), tidy = FALSE)
        }
      )

      # replacing args in the args_loop list for the subsequent matching
      args_loop[["x"]] <- prematched
    }

    # performing the matching ==================================================
    tryCatch(
      {
        verbosely(
          withr::with_preserve_seed(
            matched <- do.call(
              .match_methods[[method]]$matching_fun,
              args_loop
            )
          ),
          verbose = verbose_output
        )
      },
      error = function(e) {
        chk::abort_chk(strwrap(sprintf(
          "There was a problem with matching the samples using
        `%s`.\n Error message:
        (from `%s`) %s",
          deparse(substitute(.match_methods[[method]]$matching_fun)),
          deparse(substitute(.match_methods[[method]]$matching_fun)),
          conditionMessage(e)
        ), prefix = " ", initial = ""), tidy = FALSE)
      }
    )

    ## when no matches found, return an error
    .chk_cond(
      all(is.na(matched)),
      "No matches found!"
    )

    ## ids of the matched samples
    if (method == "nnm") {
      ids_filtered <- as.numeric(rownames(args_loop[[data_name]]))


      ids_matched <- data.frame(
        control = ids_filtered[matched$index.treated],
        treatment = ids_filtered[matched$index.control]
      )

      colnames(ids_matched) <- paste0("level_", combos[i, 1:2])
    } else if (method %in% c("fullopt", "pairopt")) {
      ids_matched <- stats::na.omit(matched)
    }

    match_results <- append(
      match_results,
      list(ids_matched)
    )
  }

  ############
  if (method == "nnm") {
    # defining the name of control column
    control_name <- paste0("level_", reference)

    # extracting the control columns
    matched_ids <- lapply(match_results, function(x) {
      sort(unique(x[, control_name, drop = TRUE]))
    })

    common_controls <- Reduce(intersect, matched_ids)

    # extract the respective matches
    extracted_matches <- lapply(match_results, function(x) {
      filter_vector <- x[, control_name] %in% common_controls
      unique(x[filter_vector, 2, drop = TRUE])
    })

    extracted_matches <- append(extracted_matches, list(common_controls))

    # merging all extracted ids
    all_extracted_ids <- Reduce(c, extracted_matches)
  } else if (method %in% c("fullopt", "pairopt")) {
    # extract only the ids of matched samples for common control identification
    matched_ids <- lapply(match_results, function(x) as.numeric(names(x)))

    # define the numbers of controls in logit_matrix
    control_range <- which(args[[data_name]]$treatment == reference)

    # extract the common controls
    matched_ids <- append(matched_ids, list(control_range))
    common_controls <- Reduce(intersect, matched_ids)

    # extract matched groups of controls
    common_controls_groups <- lapply(
      match_results,
      function(x) {
        filter_ids <- names(x) %in% common_controls
        x[filter_ids]
      }
    )

    # remove controls from match_results
    match_results <- lapply(
      match_results,
      function(x) {
        filter_ids <- names(x) %in% control_range
        x[!filter_ids]
      }
    )

    # extract matches with same groups
    extracted_matches <- mapply(function(x, y) y[y %in% x],
      common_controls_groups,
      match_results,
      SIMPLIFY = FALSE
    )

    extracted_matches <- lapply(
      extracted_matches,
      function(x) as.numeric(names(x))
    )

    extracted_matches <- append(extracted_matches, list(common_controls))

    # merging all extracted ids
    all_extracted_ids <- Reduce(c, extracted_matches)
  }

  # returning logical vector
  matched_filter <- seq_len(nrow(csmatrix)) %in% all_extracted_ids

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
