#' Title
#'
#' @param formula
#' @param data
#' @param caliper
#' @param reference
#' @param combos
#' @param kmeans.args
#' @param ratio
#' @param replace
#' @param ...
#'
#' @return
#'
#' @examples
#' @export
match_gps <- function(csmatrix,
                      caliper = 0.2,
                      reference = NULL, # not useful?
                      combos = NULL,
                      ratio = 1,
                      replace = FALSE,
                      kmeans.args = NULL, # add and process combos argument
                      kmeans.cluster = 5,
                      verbose.output = FALSE, # checking
                      ...) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  args <- list(...)
  if(is.null(kmeans.args)) kmeans.args <- list()

  # check and process the csmatrix
  if (!is.null(csmatrix)) {
    .check_df(csmatrix)
    if ("gps" %nin% class(csmatrix)) {
      "The argument `csmatrix` has to be of class `gps`."
    }

    .check_gps_matrix(csmatrix)
  } else {
    chk::abort_chk("The argument `csmatrix` is missing with no default.")
  }

  # Perform the logit transformation and combine with treatment
  logit_matrix <- cbind(treatment = csmatrix[, 1], logit(csmatrix[, -1]))

  # reference
  csmatrix$treatment <- as.factor(csmatrix$treatment)
  ref.list <- .process_ref(csmatrix$treatment,
    ordinal.treat = NULL,
    reference = reference
  )

  args[["treatment"]] <- ref.list[["data.relevel"]]
  reference <- ref.list[["reference"]]

  # kmeans.args, check list, process later
  if (!is.null(kmeans.args)) {
    chk::chk_list(kmeans.args)
  }

  # process combos
  if (is.null(combos)) {
    # generate all possible pairwise matches
    combos <- expand.grid(group1 = reference, group2 = colnames(csmatrix)[-1])
    combos <- combos[combos[, 2] != reference, ]
    #combos <- t(utils::combn(colnames(csmatrix)[-1], 2))


  } else {
    # check if data frame
    .check_df(combos, data_name = "combos")

    # check if two cols
    if (ncol(combos) != 2) {
      chk::abort_chk("The `combos` data frame must have exactly 2 columns.")
    }

    # check if all values in colnames(csmatrix)
    vectorized_combos <- as.character(c(combos[, 1], combos[, 2]))
    if (any(vectorized_combos %nin% colnames(csmatrix)[-1])) {
      chk::abort_chk("All values in the `combos` table must match the names of the unique levels of the treatment variable or the column names of the `csmatrix`.")
    }

    # check if no repeats (e.g. 1, 1)
    combos_equal <- combos[, 1] == combos[, 2]
    if (any(combos_equal)) {
      chk::abort_chk(sprintf("You can not match a group to itself, rows: %s", word_list(which(combos_equal))))
    }

    # check if all combinations unique
    # Sort each row and convert it into a string for comparison
    row_sorted <- apply(combos, 1, function(row) paste(sort(row), collapse = ","))

    # Check for duplicates
    duplicates <- duplicated(row_sorted)

    if (any(duplicates)) {
      chk::abort_chk(sprintf(
        "You can not check the same combination twice, rows: %s",
        word_list(which(duplicates))
      ))
    }
  }

  # define combos
  colnames(combos) <- c("group1", "group2")
  combos <- as.data.frame(apply(combos, 2, as.character))

  matches_n <- nrow(combos)

  # ratio
  if (length(ratio) != 1) {
    if (!.check_vecl(ratio, matches_n, check_numeric = FALSE)) {
      chk::abort_chk("The `ratio` argument must be either a single integer or an atomic vector with a length equal to the number of rows in the `combos` data frame.")
    }
  }

  .check_integer(ratio, x_name = "ratio")

  args[["M"]] <- .vectorize(ratio, matches_n)

  # replace
  if (length(replace) != 1) {
    if (!.check_vecl(replace, matches_n, check_numeric = FALSE)) {
      chk::abort_chk("The `replace` argument must be either a single logical flag or an atomic vector with a length equal to the number of rows in the `combos` data frame.")
    }

    if (!is.logical(replace) || anyNA(replace)) {
      chk::abort_chk("All values in the `replace` argument have to be logical flags.")
    }
  } else {
    chk::chk_flag(replace)
  }

  args[["replace"]] <- .vectorize(replace, matches_n)

  # process caliper
  if (is.null(caliper)) {
    chk::abort_chk("The `caliper` argument can not be NULL.")
  } else {
    if (length(caliper) != 1) {
      if (!.check_vecl(caliper, leng = matches_n, check_numeric = TRUE)) {
        chk::abort_chk("The `caliper` argument must be either a single numeric or an atomic vector with a length equal to the number of rows in the `combos` data frame.")
      }
    } else {
      chk::chk_numeric(caliper)
    }

    if (any(caliper <= 0)) {
      chk::abort_chk("The `caliper` argument has to be a positive number.")
    }
  }

  caliper <- caliper * sd(as.matrix(logit_matrix[, -1]))

  args[["caliper"]] <- .vectorize(caliper, matches_n)

  # process kmeans.cluster
  if (is.null(kmeans.cluster)) {
    chk::abort_chk("The `kmeans.cluster` argument can not be NULL.")
  } else {
    if (length(kmeans.cluster) != 1) {
      if (!.check_vecl(kmeans.cluster, leng = matches_n, check_numeric = TRUE)) {
        chk::abort_chk("The `kmeans.cluster` argument must be either a single integer or an atomic vector with a length equal to the number of rows in the `combos` data frame.")
      }
    }
  }

  .check_integer(kmeans.cluster, x_name = "kmeans.cluster")

  if (any(kmeans.cluster <= 1)) {
    chk::abort_chk("The `kmeans.cluster` argument must be an integer greater than 1.")
  }

  kmeans.args[["centers"]] <- .vectorize(kmeans.cluster, matches_n)

  ## deal with algorithm argument
  if(is.null(kmeans.args[['algorithm']])) kmeans.args['algorithm'] <- 'Hartigan-Wong'

  ## processing the kmeans arglist
  kmeans.args <- match_add_args(
    arglist = kmeans.args,
    funlist = stats::kmeans
  )

  ## vectorize the args
  kmeans.args <- lapply(kmeans.args, .vectorize, matches_n)

  ## check if tolerance specified
  tolerance_spec <- is.null(args[['tolerance']])

  ## processing the matching arglist
  args <- match_add_args(
     arglis = args,
     funlist = Matching::Matchby
  ) ## add Matching::Match()

  ## vectorize args
  args <- lapply(args, .vectorize, matches_n)
  args <- args[unlist(lapply(args, function(x) !is.null(x)))]

  ## remove tolerance if not specified earlier
  if(tolerance_spec) args <- args[names(args) %nin% 'tolerance']

  ######################## MATCHING ############################################
  ## Steps
  ## 1.) combos --> check and generate
  ## 2.) expand ratio and replace to vectors
  ## 3.) k.cluster argument --> to specify the number of clusters
  ## 4.) actual matching process
  ## 5.) merging the datasets
  ## 6.) output with special class (vector_matched)
  ## 7.) vecotrize kmeans.args and ... args
  ## 8.) refine logic for checking args

  match_results <- list()
  logit_matrix<- cbind(ID = 1:nrow(logit_matrix), logit_matrix)
  for (i in seq_len(matches_n)) {
    ## observation filter
    obs_filter <- logit_matrix[, 'treatment'] %in% as.vector(combos[i, ])

    # selecting elements from list
    kmeans.args.loop <- lapply(kmeans.args, function(x) x[[i]])
    args.loop <- lapply(args, function(x) x[[i]])

    # selecting columns for kmeans
    cols_kmeans <- colnames(logit_matrix)[-c(1:2)]
    cols_kmeans <- cols_kmeans[cols_kmeans %nin% as.vector(combos[i, ])]

    # selecting from df
    kmeans.args.loop[['x']] <- logit_matrix[, cols_kmeans]

    # fitting kmeans clusters
    tryCatch({
      verbosely(
        withr::with_preserve_seed(
          k.res <- do.call(stats::kmeans,
                           kmeans.args.loop)
        ), verbose = verbose.output
      )
    }, error = function(e) {
      chk::abort_chk(sprintf(
        "There was a problem fitting the kmeans clustering with `stats::kmeans()`.\n
                               Error message: (from `stats::kmeans()`) %s",
       conditionMessage(e)
      ), tidy = FALSE)
    })

    # filtering the clusters
    args.loop[['by']] <- k.res$cluster[obs_filter]

    # selecting columns for kmeans
    cols_matching <- colnames(logit_matrix) == combos[i, 1]

    # selecting X
    args.loop[['X']] <- logit_matrix[obs_filter, cols_matching]
    ids_df <- logit_matrix[obs_filter, 'ID']

    # defining Tr args for matching function
    args.loop[['Tr']] <- ifelse(logit_matrix[obs_filter, 2] == combos[i, 1], 0, 1)

    # performing the matching
    tryCatch({
      verbosely(
        withr::with_preserve_seed(
          matched <- do.call(Matching::Matchby,
                             args.loop)
        ), verbose = verbose.output
      )
    }, error = function(e) {
      chk::abort_chk(sprintf(
        "There was a problem with matching the samples using `Matching::Matchby()`.\n
                               Error message: (from `Matching::Matchby()`) %s",
        conditionMessage(e)
      ), tidy = FALSE)
    })

    # saving the results
    res_df <- data.frame(control = ids_df[matched$index.control],
                         treatment = ids_df[matched$index.treated])
    colnames(res_df) <- paste0('level_', combos[i, 1:2])

    match_results <- append(match_results, list(res_df))

  }

  names(match_results) <- paste0('combos_', seq_len(matches_n))

  # saving the results
  ## TO REWRITE
  ## Do this for all levels of treatment var
  ## 1.) Find all columns for given level in the match_results
  ## 2.) Pick only those ids, which are present in all results
  ## 3.) Save them
  ## 4.) Repeat for other group

  ## WORKS ONLY FOR DEFAULT COMBOS
  ## 1.) Find reference in all data frames
  ## 2.) Reduce and intersect all references
  ## 3.) filter all dataframes only for references

  # extract only refernce levels
  ref_name <- paste0('level_', reference)
  IDs_ref <- lapply(match_results, function(x) {
    x[, ref_name]
  })

  # intersection between all reference IDs
  IDs_ref_intersect <- Reduce(intersect, IDs_ref)

  # filter the list
  match_results_intersect <- lapply(match_results, function(x) {
    filter_vec <- x[, ref_name] %in% IDs_ref_intersect
    x[filter_vec, ]
  })

  # merged all results to a single vector of ids
  merged_ids <- Reduce(function(x, y) merge(x, y, by = ref_name, all = TRUE), match_results_intersect)
  merged_ids <- unname(unlist(merged_ids))

  # returning logical vector
  matched_filter <- 1:nrow(csmatrix) %in% merged_ids
  return(matched_filter)
}
