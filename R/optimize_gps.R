#' @title Optimize the Matching Process via Random Search
#'
#' @description The `optimize_gps()` function performs a random search to
#'   identify optimal combinations of parameters for the `match_gps()` and
#'   `estimate_gps()` functions. The goal is to maximize the percentage of
#'   matched samples (`perc_matched`) while minimizing the maximum standardized
#'   mean difference (`smd`), thereby improving the overall balance of
#'   covariates across treatment groups. The function supports parallel
#'   execution through the `foreach` and `future` packages, enabling
#'   multithreaded computation to accelerate the optimization process,
#'   particularly when dealing with large datasets or complex parameter spaces.
#'
#' @param data A `data.frame` containing all variables specified in the
#'   `formula` argument. If `opt_args` is used, the `data` provided within
#'   `opt_args` must match this input exactly.
#'
#' @param formula A valid formula object used to estimate the generalized
#'   propensity scores (GPS). The treatment variable appears on the left-hand
#'   side, and covariates on the right-hand side. Interactions can be specified
#'   using `*`. See [stats::formula()] and [estimate_gps()] for more details. If
#'   `opt_args` is provided, the formula within it must be identical to this
#'   argument.
#'
#' @param ordinal_treat An atomic vector defining the ordered levels of the
#'   treatment variable. This confirms the variable is ordinal and adjusts its
#'   levels accordingly using
#'   `factor(treat, levels = ordinal_treat, ordered = TRUE)`. It is passed
#'   directly to `estimate_gps()`. If `NULL`, ordinal GPS estimation methods
#'   such as `polr` will be excluded from the optimization. See [estimate_gps()]
#'   for details.
#'
#' @param n_iter Integer. Number of unique parameter combinations to evaluate
#'   during optimization. Higher values generally yield better results but
#'   increase computation time. For large datasets or high-dimensional parameter
#'   spaces, increasing `n_iter` is recommended. When using parallel processing
#'   (`n_cores > 1`), performance gains become more apparent with larger
#'   `n_iter`. Too many cores with too few iterations may introduce
#'   overhead and reduce efficiency.
#'
#' @param opt_args An object of class `"opt_args"` containing optimization
#'   parameters and argument settings. Use [make_opt_args()] to create this
#'   object. It specifies the search space for the GPS estimation and matching
#'   procedure.
#'
#' @details The output is an S3 object of class `best_opt_result`. Its core
#'   component is a `data.frame` containing the parameter settings for the
#'   best-performing models, grouped and ranked based on their balance quality.
#'
#' Optimization results are categorized into seven bins based on the maximum
#' standardized mean difference (SMD):
#' - 0.00-0.05
#' - 0.05-0.10
#' - 0.10-0.15
#' - 0.15-0.20
#' - 0.20-0.25
#' - 0.25-0.30
#' - Greater than 0.30
#'
#' Within each SMD group, the parameter combination(s) achieving the highest
#' `perc_matched` (i.e., percentage of matched samples) is selected. In cases
#' where multiple combinations yield identical `smd` and `perc_matched`, all
#' such results are retained. Combinations where matching failed or GPS
#' estimation did not converge will return `NA` in the result columns (e.g.,
#' `perc_matched`, `smd`).
#'
#' The returned `data.frame` includes the following columns (depending on the
#' number of treatment levels):
#' - `iter_ID`: Unique identifier for each parameter combination
#' - `method_match`: Matching method used in [match_gps()], e.g., `"nnm"` or
#' `"fullopt"`
#' - `caliper`: Caliper value used in [match_gps()]
#' - `order`: Ordering of GPS scores prior to matching
#' - `kmeans_cluster`: Number of k-means clusters used
#' - `replace`: Whether replacement was used in matching (`nnm` only)
#' - `ties`: Tie-breaking rule in nearest-neighbor matching (`nnm` only)
#' - `ratio`: Control-to-treated ratio for `nnm`
#' - `min_controls`, `max_controls`: Minimum and maximum controls for `fullopt`
#' - `reference`: Reference group used in both [estimate_gps()] and
#' [match_gps()]
#' - `perc_matched`: Percentage of matched samples (from [balqual()])
#' - `smd`: Maximum standardized mean difference (from [balqual()])
#' - `p_{group_name}`: Percent matched per treatment group (based on group s
#' ample size)
#' - `method_gps`: GPS estimation method used (from [estimate_gps()])
#' - `link`: Link function used in GPS model
#' - `smd_group`: SMD range category for the row
#'
#' The resulting `best_opt_result` object also includes a custom `summary()`
#' method that summarizes:
#' - The number of optimal parameter sets per SMD group
#' - Their associated SMD and match rates
#' - Total combinations tested
#' - Total runtime of the optimization loop
#'
#'
#' @return An S3 object of class `best_opt_result`. The core component is a
#'   `data.frame` containing the parameter combinations and results of the
#'   optimization procedure. You can access it using `attr(result,
#'   "opt_results")` or by calling `View(result)`, where `result` is your
#'   `best_opt_result` object.
#'
#' The object contains the following custom attributes:
#'
#' - **`opt_results`**: A `data.frame` of optimization results. Each row
#' corresponds to a unique parameter combination tested. For a complete
#' description of columns, see the ***Details*** section.
#'
#' - **`optimization_time`**: Time (in seconds) taken by the optimization loop
#' (i.e., the core `for`-loop that evaluates combinations). This does **not**
#' include the time needed for GPS estimation, pre-processing, or merging of
#' results after loop completion. On large datasets, these excluded steps can
#' still be substantial.
#'
#' - **`combinations_tested`**: Total number of unique parameter combinations
#' evaluated during optimization.
#'
#' - **`smd_results`**: A detailed table of standardized mean differences (SMDs)
#'  for all pairwise treatment group comparisons and for all covariates
#'  specified in the `formula`. This is used by the [select_opt()] function to
#'  filter optimal models based on covariate-level balance across groups.
#'
#' - **`treat_names`**: A character vector with the names of the unique
#' treatment groups.
#'
#' - **`model_covs`**: A character vector listing the model covariates (main
#' effects and interactions) used in the `formula`. These names correspond to
#' the variables shown in the `smd_results` table.
#'
#' @examples
#' # Define formula for GPS estimation and matching
#' formula_cancer <- formula(status ~ age * sex)
#'
#' # Set up the optimization parameter space
#' opt_args <- make_opt_args(cancer, formula_cancer, gps_method = "m1")
#'
#' # Run optimization with 2000 random parameter sets and a fixed seed
#' \donttest{
#' withr::with_seed(
#'   8252,
#'   {
#'     optimize_gps(
#'       data = cancer,
#'       formula = formula_cancer,
#'       opt_args = opt_args,
#'       n_iter = 2000
#'     )
#'   }
#' )
#' }
#' @export
optimize_gps <- function(data = NULL,
                         formula,
                         ordinal_treat = NULL,
                         n_iter = 1000,
                         opt_args = NULL) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################

  # check that the n_iter is integer and greater than 1
  .check_integer(n_iter)
  .chk_cond(
    n_iter < 1,
    "The number of iterations (n_iter) has to be at least 1."
  )

  # -------- parallel backend ---------------------------------------------
  # NOTE: the parallel backend for foreach (e.g. via the future or
  # doFuture ecosystems) is expected to be registered *outside* of
  # this function. For example:
  #
  # future::plan(future::multisession, workers = 5)
  # doFuture::registerDoFuture()
  # doRNG::registerDoRNG(once = FALSE)
  #
  # Here we only:
  # (i) ensure required packages for the loops are available and
  # (ii) choose between sequential (%do%) and parallel (%dorng%)
  #     execution based on the currently registered foreach backend.

  # packages required for the foreach + progressr loops
  rlang::check_installed(
    c("foreach", "progressr"),
    reason = "to execute optimization loops"
  )

  # detect number of workers from the registered foreach backend
  n_workers <- foreach::getDoParWorkers()
  if (is.null(n_workers) || !is.finite(n_workers)) {
    n_workers <- 1L
  }

  if (n_workers <= 1L) {
    # sequential - only foreach needed
    rlang::check_installed(
      "foreach",
      reason = "for sequential foreach execution"
    )

    # get the '%do%' operator from the foreach namespace
    `%doparallel%` <- get("%do%", asNamespace("foreach"))

  } else {
    # parallel - foreach + doRNG needed
    rlang::check_installed(
      c("doRNG"),
      reason = "for reproducible parallel foreach execution"
    )

    # get the '%dorng%' operator from the doRNG namespace
    `%doparallel%` <- get("%dorng%", asNamespace("doRNG"))
  }

  ## adding default optimization parameters if not specified
  if (is.null(opt_args)) {
    opt_args <- make_opt_args(
      data    = data,
      formula = formula
    )
  }

  # validate class of opt_args if provided
  .chk_cond(
    !inherits(opt_args, "opt_args"),
    "opt_args must be of class 'opt_args'. Use make_opt_args() to create it."
  )

  ###################### DEFINING SEARCH SPACE #################################
  ###### ESTIMATE GPS

  ## define the estimate_space for the function estimate_gps

  # defining available gps methods from the .gps_methods object
  available_methods <- names(.gps_methods)

  # create method-link combinations using lapply + do.call
  # (multiple links for one method available)
  estimate_methods <- do.call(
    rbind,
    lapply(
      available_methods,
      function(method) {
        links <- .gps_methods[[method]]$link_fun
        data.frame(method = method, link = links)
      }
    )
  )

  # define the short name for the method
  estimate_methods$gps_method <- paste0("m", 1:10) ## --> see make_opt_args()!

  # subset for gps_method specified in the opt_args
  estimate_methods <- estimate_methods[
    estimate_methods$gps_method %in% opt_args[["gps_method"]],
  ]

  # Create reference levels from unique treatment values
  available_refs <- data.frame(refs = opt_args[["reference"]])

  # Cartesian product of methods x reference levels
  estimate_space <- merge(
    estimate_methods,
    available_refs,
    by = NULL
  )

  # Adding unique names to each column
  estimate_space$row_name <- paste(
    rep("estimate", nrow(estimate_space)),
    seq_len(nrow(estimate_space)),
    sep = "_"
  )

  # remove "polr" from the estimate space if ordinal_treat is NULL
  if (is.null(ordinal_treat)) {
    estimate_space <- estimate_space[estimate_space$method != "polr", ]
  }

  ##############################################################################
  ###### MATCH GPS

  ## defining number of iterations and searches
  n_iter_final <- n_iter # final number of iterations used for looping
  n_iter <- n_iter * 10 # number of iterations for grid search

  # generate base parameter grid using random sampling
  withr::with_preserve_seed({
    search_matching <- data.frame(
      gps_model = sample(
        unique(estimate_space$row_name),
        n_iter,
        replace = TRUE
      ),
      method = sample(
        opt_args[["matching_method"]],
        n_iter,
        replace = TRUE
      ),
      caliper = sample(
        opt_args[["caliper"]],
        n_iter,
        replace = TRUE
      ),
      order = sample(
        opt_args[["order"]],
        n_iter,
        replace = TRUE
      ),
      kmeans_cluster = sample(
        opt_args[["cluster"]],
        n_iter,
        replace = TRUE
      )
    )
  })

  # Preallocate new columns to avoid growing the dataframe in a loop and
  # dimension mismatch
  cols_to_add <- c("replace", "ties", "ratio", "min_controls", "max_controls")
  search_matching[cols_to_add] <- lapply(cols_to_add, function(x) NA)

  # Logical index for method-specific assignments
  is_nnm <- search_matching$method == "nnm"
  is_fullopt <- search_matching$method == "fullopt"

  # Assign values for "nnm" method
  n_nnm <- sum(is_nnm)
  withr::with_preserve_seed({
    search_matching$replace[is_nnm] <- sample(
      opt_args[["replace"]],
      n_nnm,
      replace = TRUE
    )
  })
  withr::with_preserve_seed({
    search_matching$ties[is_nnm] <- sample(
      opt_args[["ties"]],
      n_nnm,
      replace = TRUE
    )
  })
  withr::with_preserve_seed({
    search_matching$ratio[is_nnm] <- sample(
      opt_args[["ratio"]],
      n_nnm,
      replace = TRUE
    )
  })

  # Assign values for "fullopt"
  n_fullopt <- sum(is_fullopt)
  withr::with_preserve_seed({
    search_matching$min_controls[is_fullopt] <- sample(
      opt_args[["min_controls"]],
      n_fullopt,
      replace = TRUE
    )
  })
  withr::with_preserve_seed({
    search_matching$max_controls[is_fullopt] <- sample(
      opt_args[["max_controls"]],
      n_fullopt,
      replace = TRUE
    )
  })

  # note: all sample() calls are wrapped inside with_preserve_seed --> without
  # it they somehow managed to change the global user seed

  # add reference from estimate_space
  search_matching <- merge(
    search_matching,
    estimate_space[, c("refs", "row_name")],
    by.x = "gps_model",
    by.y = "row_name",
    all.x = TRUE
  )

  # Changing refs to reference
  names(search_matching)[names(search_matching) == "refs"] <- "reference"
  search_matching$reference <- as.character(search_matching$reference)

  # quality control

  ## remove duplicates
  search_matching <- search_matching[!duplicated(search_matching), ]

  ## ensure that max_controls >= min_controls
  if ("fullopt" %in% opt_args[["matching_method"]]) {
    valid <- search_matching$min_controls <= search_matching$max_controls
    valid[is.na(valid)] <- TRUE
    search_matching <- search_matching[valid, ]
  }

  ## pick final number of samples
  search_matching <- tryCatch(
    {
      withr::with_preserve_seed({
        search_matching[
          sample(
            nrow(search_matching),
            n_iter_final,
            replace = FALSE
          ),
        ]
      })
    },
    error = function(e) {
      warning(
        paste(
          "Sampling the final number of rows without replacement failed.",
          "Consider expanding the parameter range.",
          "Falling back to sampling with replacement."
        )
      )
      withr::with_preserve_seed({
        search_matching[
          sample(
            nrow(search_matching),
            n_iter_final,
            replace = TRUE
          ),
        ]
      })
    }
  )

  ###################### FITTING ESTIMATE GPS #################################

  # CLI-style message
  cli::cli_alert_info("Initiating estimation of the GPS...")

  ## the problem with seeds in the optimization function is, that we
  ## need EACH ITERATION to run with the EXACT SAME SEED. I think the
  ## easiest way to achieve this is save the current .Random.seed,
  ## then export it and use withr::with_seed() to make sure, that
  ## the seed remains exactly the same for each i

  export_seed <- .Random.seed

  with_rng_state <- function(seed_state, code) {
    old_state <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }

    assign(".Random.seed", seed_state, envir = .GlobalEnv)

    on.exit(
      {
        if (is.null(old_state)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else {
          assign(".Random.seed", old_state, envir = .GlobalEnv)
        }
      },
      add = TRUE
    )

    force(code)
  }

  # Enable global handler for the progress (once is enough)
  progressr::handlers("txtprogressbar")

  # to avoid seed leaks
  withr::with_preserve_seed({
    ## defining the loop length
    loop_seq <- seq_len(nrow(estimate_space))

    ## looping through all estimate space combinations
    suppressMessages({
      estimate_results <- foreach::foreach(
        i = loop_seq,
        .packages = c("vecmatch", "withr"), # add any other needed packages
        # i needed to export all used objects to the workers in the parallel
        # setup
        .export = c(
          "data", "formula", "estimate_space", "ordinal_treat",
          "estimate_gps", "csregion", "export_seed"
        ),
        .errorhandling = "pass"
      ) %doparallel% {
        # i know it may seems unnecessary to wrap it all inside a function and
        # call it at the end, but it has a purpose. During the debugging i
        # noticed that the workers somehow share a common environment and
        # managed to change the values of some variables during the run. It
        # was a total mess.
        #
        # so the solution was to wrap the whole loop inside a function and
        # call it at the end. Firstly, it creates an isolated environment for
        # each iteration (function env), and the whole iteration is executed
        # only within this environment, which is then deleted

        run_iteration <- function(i) {
          # defining the current argument list
          arglist_loop <- list(
            data          = data,
            formula       = formula,
            method        = estimate_space[i, "method"],
            link          = estimate_space[i, "link"],
            reference     = as.character(estimate_space[i, "refs"]),
            ordinal_treat = ordinal_treat
          )

          # estimating the gps
          # double caution - better safe than sorry
          withr::with_preserve_seed({
            loop_estimate <- tryCatch(
              {
                do.call(estimate_gps, arglist_loop)
              },
              error = function(e) {
                list(
                  error = TRUE,
                  message = sprintf(
                    "Error with method %s: %s",
                    arglist_loop$method,
                    conditionMessage(e)
                  )
                )
              }
            )

            # calculating csregion borders
            loop_estimate <- csregion(loop_estimate)
          })

          # remove objects
          rm(list = setdiff(ls(), "loop_estimate"))

          # defining the output
          loop_estimate <- list(loop_estimate)
          names(loop_estimate) <- estimate_space[i, "row_name"]

          loop_estimate
        }

        # call the function within isolated env (function env)
        with_rng_state(export_seed, {
          run_iteration(i)
        })
      }
    })
  })

  ## unnesting the estimate_results list
  list_names <- vapply(
    estimate_results,
    function(x) names(x[1]),
    character(1)
  )
  estimate_results <- lapply(estimate_results, function(x) x[[1]])

  # assigning unique row_names for the estimates
  names(estimate_results) <- list_names

  ###################### MATCHING OPTIMIZER ####################################
  # The search grid is already generated

  # CLI-style messages
  cli::cli_alert_info("Initiating matching optimization...")
  cli::cli_alert_info("Optimizing matching...")

  ## starting time
  time_start <- Sys.time()

  ## preserving the seed (avoiding seed leaks)
  withr::with_preserve_seed({
    suppressMessages({
      progressr::with_progress({
        ## loop length and progress bar
        loop_seq <- seq_len(n_iter_final)
        throttle <- 100L

        ## number of progress steps:
        ## - ceiling() so we NEVER have 0 steps
        ## - max(1L, ...) as extra safety
        n_steps <- max(1L, ceiling(n_iter_final / throttle))

        ## one progressor for this loop
        p <- progressr::progressor(steps = n_steps)

        # Precompute unique treatment names outside the foreach loop
        all_treatments <- unique(estimate_results[[1]]$treatment)
        treatment_cols <- paste0("p_", all_treatments)

        smd_colnames <- c("group1", "group2", attr(opt_args, "model_covs"))
        smd_template <- as.data.frame(
          matrix(NA, nrow = 1, ncol = length(smd_colnames)),
          stringsAsFactors = FALSE
        )
        colnames(smd_template) <- smd_colnames

        opt_results <- foreach::foreach(
          i = loop_seq,
          .packages = "vecmatch",
          .export = c(
            "search_matching", "estimate_results", "formula",
            "treatment_cols", "smd_colnames", "smd_template",
            "match_gps", "balqual", "export_seed", # make sure workers see these
            "p", "throttle", "loop_seq"
          ),
          .errorhandling = "pass"
        ) %doparallel% {
          run_iteration <- function(i) {
            ## define iter ID
            iter_ID <- paste0("ID", i)

            ## processing the argument list
            args_loop <- as.list(search_matching[i, , drop = FALSE])
            args_loop <- args_loop[!is.na(args_loop)]
            args_loop[["csmatrix"]] <- estimate_results[[args_loop$gps_model]]
            args_loop <- args_loop[-1] # drop gps_model

            # Defaults in case of error
            perc_matched <- NA_real_
            smd <- NA_real_
            perc <- as.data.frame(
              t(rep(NA_real_, length(treatment_cols)))
            )
            colnames(perc) <- treatment_cols

            smd_df <- smd_template # predefine smd_df with correct structure

            try(
              {
                # matching
                loop_estimate <- do.call(match_gps, args_loop)

                # max SMD and %matched statistics
                utils::capture.output({
                  qual_out <- balqual(
                    loop_estimate,
                    formula,
                    type      = "smd",
                    statistic = "max",
                    round     = 8,
                    print_out = FALSE
                  )
                })

                # Take the smd_df from balqual
                smd_extracted <- attr(qual_out, "smd_df_combo")

                ## baseline stats
                perc_matched <- qual_out$perc_matched
                smd <- qual_out$summary_head

                # Calculate percentages
                ptab <- as.data.frame(qual_out$count_table)
                ptab$Before <- as.numeric(ptab$Before)
                ptab$After <- as.numeric(ptab$After)
                ptab$p <- (ptab$After / ptab$Before) * 100

                # Fill into correct named columns
                computed <- stats::setNames(
                  as.list(ptab$p),
                  paste0("p_", ptab$Treatment)
                )
                for (col in names(computed)) {
                  if (col %in% treatment_cols) {
                    perc[[col]] <- computed[[col]]
                  }
                }

                # update template with correct number of rows
                smd_df <- as.data.frame(
                  matrix(
                    NA,
                    nrow = nrow(smd_extracted),
                    ncol = length(smd_colnames)
                  ),
                  stringsAsFactors = FALSE
                )
                colnames(smd_df) <- smd_colnames

                # fill only matching columns
                for (col in colnames(smd_extracted)) {
                  if (col %in% smd_colnames) {
                    smd_df[[col]] <- smd_extracted[[col]]
                  }
                }
              },
              silent = TRUE
            )

            ## Update progress at most n_steps times
            if (i %% throttle == 0L) {
              p(sprintf("Running %d/%d", i, max(loop_seq)))
            }

            # setting up the resulting data frame
            result_row <- cbind(
              iter_ID       = iter_ID,
              search_matching[i, ],
              perc_matched  = perc_matched,
              smd           = smd,
              perc
            )

            # corresponding smd data.frame
            smd_df <- cbind(iter_ID = iter_ID, smd_df)

            # returning output from single iteration
            res <- list(
              results = as.data.frame(result_row, stringsAsFactors = FALSE),
              smd_dfs = as.data.frame(smd_df, stringsAsFactors = FALSE)
            )

            ## --- tidy-up: remove *everything* except 'res' ----------------
            rm(list = setdiff(ls(), "res"))
            ## optional: a *light* sweep every 256th iteration
            # if (i %% 1024 == 0L) gc(FALSE)

            # returning the results
            res
          }

          # running the function inside an isolated env
          with_rng_state(export_seed, {
            run_iteration(i)
          })
        }
      })
    })
  })

  # defining the stop time and running time
  time_stop <- Sys.time()
  time_diff <- round(as.numeric(time_stop - time_start, units = "secs"), 2)

  ## extract looping results
  opt_results_df <- do.call(
    rbind,
    lapply(opt_results, `[[`, "results")
  )
  smd_df_all <- do.call(
    rbind,
    lapply(opt_results, `[[`, "smd_dfs")
  )

  # processing the results

  ## merge the gps_model
  opt_results_df <- merge(
    opt_results_df,
    estimate_space[, c("method", "link", "row_name")],
    by.x = "gps_model",
    by.y = "row_name"
  )

  ## rename method.y and method.x
  names(opt_results_df)[names(opt_results_df) == "method.y"] <- "method_gps"
  names(opt_results_df)[names(opt_results_df) == "method.x"] <- "method_match"

  # Define SMD intervals
  breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, Inf)
  labels <- c(
    "0-0.05", "0.05-0.10", "0.10-0.15",
    "0.15-0.20", "0.20-0.25", "0.25-0.30",
    ">0.30"
  )

  # Cut SMD into intervals
  opt_results_df$smd_group <- cut(
    opt_results_df$smd,
    breaks  = breaks,
    labels  = labels,
    right   = FALSE
  )

  ### BEST RESULTS OUTPUT FILTERING SECTION ####################################

  # Remove rows with missing perc_matched or smd_group
  filtered <- opt_results_df[
    !is.na(opt_results_df$perc_matched) &
      !is.na(opt_results_df$smd_group),
  ]

  # Get unique SMD groups
  groups <- unique(filtered$smd_group)

  # Initialize empty list to collect best rows
  best_rows_list <- list()

  # Loop through each group and extract rows with max perc_matched
  for (g in groups) {
    group_rows <- filtered[filtered$smd_group == g, ]
    max_value <- max(group_rows$perc_matched, na.rm = TRUE)
    best_rows <- group_rows[group_rows$perc_matched == max_value, ]
    best_rows_list[[as.character(g)]] <- best_rows
  }

  # Combine all best rows into one data.frame
  best_rows_final <- do.call(rbind, best_rows_list)

  ## reset row_names and remove gps_model + ID
  rownames(best_rows_final) <- NULL
  best_rows_final <- best_rows_final[, -1]

  # assemble result object: attributes + class in ONE place
  results <- structure(
    best_rows_final,
    optimization_time   = time_diff,
    combinations_tested = format(n_iter_final, scientific = FALSE),
    opt_results         = opt_results_df,
    smd_results         = smd_df_all,
    treat_names         = unique(estimate_results[[1]]$treatment),
    model_covs          = attr(opt_args, "model_covs"),
    function_call       = match.call(),
    class               = c("best_opt_result", "data.frame")
  )

  invisible(results)
}

#' @export
summary.best_opt_result <- function(object, digits = 3, ...) {
  # keep attributes BEFORE any subsetting (subsetting a data.frame drops them)
  time_taken <- attr(object, "optimization_time", exact = TRUE)
  n_combos <- attr(object, "combinations_tested", exact = TRUE)

  x <- object

  # require necessary columns
  needed <- c("smd", "perc_matched", "smd_group")

  # filter to rows with non-NA smd / perc_matched / smd_group
  ok <- !is.na(x$smd) & !is.na(x$perc_matched) & !is.na(x$smd_group)
  x <- x[ok, , drop = FALSE]

  # split by SMD group
  groups <- split(x, x$smd_group, drop = TRUE)

  # order groups by minimal SMD (lower SMD first)
  smd_order <- vapply(
    groups,
    function(g) {
      min(g$smd, na.rm = TRUE)
    },
    numeric(1L)
  )
  sorted_groups <- names(sort(smd_order, decreasing = FALSE))

  # build summary table (same logic as before, but without printing)
  result_table <- data.frame(
    smd_group = character(),
    unique_configs = integer(),
    smd = numeric(),
    perc_matched = numeric(),
    stringsAsFactors = FALSE
  )

  for (grp in sorted_groups) {
    rows <- groups[[grp]]

    # columns that define a configuration
    config_cols <- setdiff(names(rows), c("smd", "perc_matched"))
    unique_configs <- nrow(unique(rows[, config_cols, drop = FALSE]))

    smd_vals <- rows$smd
    perc_vals <- rows$perc_matched

    summary_rows <- unique(data.frame(
      smd_group = grp,
      unique_configs = unique_configs,
      smd = smd_vals,
      perc_matched = perc_vals,
      stringsAsFactors = FALSE
    ))

    result_table <- rbind(result_table, summary_rows)
  }


  # return a proper summary object
  res <- structure(
    result_table,
    optimization_time   = time_taken,
    combinations_tested = n_combos,
    digits              = as.integer(digits),
    class               = c("summary.best_opt_result", "data.frame")
  )

  res
}

#' @export
print.summary.best_opt_result <- function(x,
                                          digits = attr(x, "digits",
                                            exact = TRUE
                                          ) %||% 3L,
                                          ...) {
  digits <- as.integer(digits)
  if (!is.finite(digits) || digits < 0L) {
    digits <- 3L
  }

  cat("Best Optimization Results by SMD Group\n")
  cat("======================================\n\n")

  result_table <- x

  # prepare table header
  header <- c("smd_group", "unique_configs", "smd", "perc_matched")
  col_widths <- c(12, 16, 8, 14)
  total_width <- sum(col_widths) + length(col_widths) + 1

  cat(strrep("-", total_width), "\n", sep = "")
  cat(
    "|",
    paste(mapply(function(name, width) {
      format(name, width = width, justify = "centre")
    }, header, col_widths), collapse = "|"),
    "|\n",
    sep = ""
  )
  cat(strrep("-", total_width), "\n", sep = "")

  for (i in seq_len(nrow(result_table))) {
    row <- result_table[i, ]
    cat(
      "|",
      format(row$smd_group, width = col_widths[1], justify = "left"), "|",
      format(row$unique_configs,
        width    = col_widths[2],
        justify  = "right"
      ), "|",
      format(round(row$smd, digits),
        width    = col_widths[3],
        justify  = "right"
      ), "|",
      format(round(row$perc_matched, digits),
        width    = col_widths[4],
        justify  = "right"
      ), "|\n",
      sep = ""
    )
  }


  cat(strrep("-", total_width), "\n\n", sep = "")

  # Optimization summary footer
  time_taken <- attr(x, "optimization_time", exact = TRUE)
  n_combos <- attr(x, "combinations_tested", exact = TRUE)

  if (!is.null(time_taken) || !is.null(n_combos)) {
    cat("Optimization Summary\n")
    cat("--------------------\n")
  }

  if (!is.null(n_combos)) {
    cat(
      "Total combinations tested  : ",
      format(n_combos, big.mark = ""),
      "\n",
      sep = ""
    )
  }

  if (!is.null(time_taken)) {
    cat(
      "Total optimization time [s]: ",
      format(time_taken),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}

# internal helper: print best_opt_result with a nice header
.print_best_opt_result_core <- function(x, n_show = 10L, ...) {
  # x: best_opt_result object (data.frame with optimization summary)

  n <- nrow(x)
  p <- ncol(x)

  opt_time <- attr(x, "optimization_time", exact = TRUE)
  comb_test <- attr(x, "combinations_tested", exact = TRUE)
  opt_res <- attr(x, "opt_results", exact = TRUE)
  smd_res <- attr(x, "smd_results", exact = TRUE)
  treat_attr <- attr(x, "treat_names", exact = TRUE)
  model_covs <- attr(x, "model_covs", exact = TRUE)

  # coerce treatments to character for printing
  if (!is.null(treat_attr)) {
    if (is.factor(treat_attr)) {
      treat_levels <- levels(treat_attr)
    } else {
      treat_levels <- unique(as.character(treat_attr))
    }
  } else {
    treat_levels <- character(0)
  }

  # basic SMD / perc summary
  smd_vals <- x[["smd"]]
  pm_vals <- x[["perc_matched"]]
  smd_range <- if (all(is.na(smd_vals))) {
    "<all NA>"
  } else {
    sprintf("[%.3f, %.3f]", min(smd_vals, na.rm = TRUE),
            max(smd_vals, na.rm = TRUE))
  }
  pm_range <- if (all(is.na(pm_vals))) {
    "<all NA>"
  } else {
    sprintf("[%.1f, %.1f]", min(pm_vals, na.rm = TRUE),
            max(pm_vals, na.rm = TRUE))
  }

  cli::cli_ul()
  cli::cli_li("Rows (selected configurations): {n}")
  cli::cli_li("Columns: {p}")
  cli::cli_li("Optimization time (sec): {opt_time %||% '<unknown>'}")
  cli::cli_li("Combinations tested: {comb_test %||% '<unknown>'}")

  cli::cli_li("Treatments: {.field {if (length(treat_levels))
              paste(treat_levels, collapse = ', ') else '<none>'}}")
  cli::cli_li("Number of covariates in balance check:
              {length(model_covs %||% character(0))}")

  cli::cli_li("SMD range in selected set: {smd_range}")
  cli::cli_li("% matched range in selected set: {pm_range}")
  cli::cli_end()

  cli::cli_text("")

  n_show <- as.integer(n_show)
  if (!is.finite(n_show) || n_show <= 0L) {
    n_show <- 10L
  }

  if (n > n_show) {
    cli::cli_text("Showing the first {n_show} of {n} rows:")
  } else {
    cli::cli_text("Showing all {n} rows:")
  }

  cli::cli_text("")

  base::print.data.frame(utils::head(x, n_show), ...)

  invisible(x)
}

#' @export
print.best_opt_result <- function(x, ...) {
  cli::cli_text("{.strong best_opt_result object} (GPS matching
                optimization summary)")
  cli::cli_text("")
  .print_best_opt_result_core(x, ...)
}

#' @export
plot.best_opt_result <- function(x,
                                 smd_xlim = c(0, 1),
                                 smd_cutoffs = c(0.10, 0.25),
                                 xlab = "Max SMD",
                                 ylab = "% matched",
                                 main = "Matching optimization results",
                                 ...) {
  # x: best_opt_result object (data.frame with columns 'smd' and 'perc_matched')

  if (!("smd" %in% names(x) && "perc_matched" %in% names(x))) {
    stop(
      "Object of class 'best_opt_result' must contain columns ",
      "'smd' and 'perc_matched'."
    )
  }

  smd <- x[["smd"]]
  pm <- x[["perc_matched"]]

  # drop NAs
  ok <- !is.na(smd) & !is.na(pm)
  smd <- smd[ok]
  pm <- pm[ok]

  if (!length(smd)) {
    warning("No finite (smd, perc_matched) pairs to plot.")
    return(invisible(x))
  }

  # enforce numeric
  smd <- as.numeric(smd)
  pm <- as.numeric(pm)

  # x-axis limited to [0, 1] as requested
  smd_xlim <- c(0, 1)

  # a reasonable default for y-axis: 0–100 (percent)
  pm_ylim <- range(c(0, 100, pm), na.rm = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::plot(
    smd,
    pm,
    xlim = smd_xlim,
    ylim = pm_ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    pch  = 19,
    ...
  )

  # vertical cutoff lines
  smd_cutoffs <- smd_cutoffs[is.finite(smd_cutoffs)]
  for (v in smd_cutoffs) {
    if (v >= smd_xlim[1] && v <= smd_xlim[2]) {
      graphics::abline(v = v, lty = 2)
    }
  }

  # optional legend for cutoffs (only if any cutoffs are inside range)
  if (length(smd_cutoffs)) {
    graphics::legend(
      "bottomleft",
      legend = paste0("SMD = ", smd_cutoffs),
      lty    = 2,
      bty    = "n"
    )
  }

  invisible(x)
}


#' @export
str.best_opt_result <- function(object, ...) {
  n <- nrow(object)
  p <- ncol(object)

  opt_time <- attr(object, "optimization_time", exact = TRUE)
  comb_test <- attr(object, "combinations_tested", exact = TRUE)
  opt_res <- attr(object, "opt_results", exact = TRUE)
  smd_res <- attr(object, "smd_results", exact = TRUE)
  treat_attr <- attr(object, "treat_names", exact = TRUE)
  model_covs <- attr(object, "model_covs", exact = TRUE)
  call <- attr(object, "function_call", exact = TRUE)

  ## Treatments -----------------------------------------------------------
  if (!is.null(treat_attr)) {
    if (is.factor(treat_attr)) {
      treat_levels <- levels(treat_attr)
    } else {
      treat_levels <- unique(as.character(treat_attr))
    }
  } else {
    treat_levels <- character(0)
  }

  k_treat <- length(treat_levels)

  ## Header --------------------------------------------------------------
  cat("best_opt_result object: GPS matching optimization summary\n")
  cat(sprintf(" Dimensions: %d rows x %d columns\n", n, p))

  cat(sprintf(
    " Optimization time (sec): %s\n",
    if (!is.null(opt_time)) as.character(opt_time) else "<unknown>"
  ))
  cat(sprintf(
    " Combinations tested: %s\n",
    if (!is.null(comb_test)) as.character(comb_test) else "<unknown>"
  ))

  cat(sprintf(
    " Number of treatment levels: %d\n",
    k_treat
  ))
  cat(sprintf(
    " Treatment levels: %s\n",
    if (k_treat) paste(treat_levels, collapse = ", ") else "<none>"
  ))

  cat(sprintf(
    " Number of covariates in balance check: %d\n",
    length(model_covs %||% character(0))
  ))
  if (!is.null(model_covs)) {
    cat(sprintf(
      " Covariates: %s\n",
      paste(model_covs, collapse = ", ")
    ))
  }

  if (!is.null(opt_res)) {
    cat(sprintf(
      " opt_results: data.frame with %d rows and %d columns\n",
      NROW(opt_res), NCOL(opt_res)
    ))
  } else {
    cat(" opt_results: <none>\n")
  }

  if (!is.null(smd_res)) {
    cat(sprintf(
      " smd_results: data.frame with %d rows and %d columns\n",
      NROW(smd_res), NCOL(smd_res)
    ))
  } else {
    cat(" smd_results: <none>\n")
  }

  if (!is.null(call)) {
    cat(" Call:\n")
    cat(
      "  ",
      paste(deparse(call, width.cutoff = 80L), collapse = "\n  "),
      "\n",
      sep = ""
    )
  }

  cat("\nUnderlying data.frame structure:\n")

  ## Delegate to data.frame method for the actual structure -------------
  utils::str(unclass(object), ...)

  invisible(object)
}

#' @title Select Optimal Parameter Combinations from Optimization Results
#'
#' @description `select_opt()` is a helper function to filter and prioritize
#' results from `optimize_gps()` based on the specific goals of a study.
#' Depending on the research design, certain pairwise comparisons or treatment
#' groups may be more important than others. For example:
#'
#' - You may want to prioritize matching between a specific groups (e.g.
#' specific disease vs. controls), while ignoring other group comparisons
#' during SMD evaluation.
#' - You may wish to retain as many samples as possible from a critical group or
#'  set of groups, regardless of matching rates in other groups.
#'
#' This function enables targeted selection of optimal parameter combinations
#' by:
#' - Evaluating SMDs for specific pairwise treatment group comparisons,
#' - Selecting key covariates to assess balance,
#' - Prioritizing matched sample size in selected treatment groups.
#'
#' By combining these criteria, `select_opt()` allows you to tailor the
#' optimization output to your study's focus - whether it emphasizes covariate
#' balance in targeted group comparisons or maximizing sample retention for
#' specific subgroups.
#'
#' @param x An object of class `best_opt_result`, produced by the
#'   `optimize_gps()` function.
#'
#' @param smd_groups A `list` of pairwise comparisons (as `character` vectors of
#'   length 2) specifying which treatment group comparisons should be
#'   prioritized in SMD evaluation. Each element must be a valid pair of
#'   treatment levels. If `NULL`, all pairwise comparisons are used. Example:
#'   `list(c("adenoma", "crc_malignant"), c("controls", "adenoma"))`
#'
#' @param smd_variables A `character` vector of covariate names to include in
#'   the SMD evaluation. Must match variables listed in `attr(x, "model_covs")`.
#'
#' @param smd_type A `character` string (`"mean"` or `"max"`), defining how to
#'   aggregate SMDs across covariates and comparisons. `"max"` selects
#'   combinations with the lowest maximum SMD; `"mean"` uses the average SMD.
#'
#' @param perc_matched A `character` vector of treatment levels for which the
#'   matching rate should be maximized. If `NULL`, overall `perc_matched` is
#'   used. If specified, only the sum of matching percentages for the listed
#'   groups is used for selection within each SMD category.
#'
#' @details Optimization results are grouped into bins based on the
#' **maximum SMD** observed for each parameter combination. These bins follow
#' the same structure as in `optimize_gps()`:
#'
#' - 0.00-0.05
#' - 0.05-0.10
#' - 0.10-0.15
#' - 0.15-0.20
#' - 0.20-0.25
#' - 0.25-0.30
#' - 0.30-0.35
#' - 0.35-0.40
#' - 0.40-0.45
#' - 0.45-0.50
#' - more than 0.50
#'
#' Within each bin, models are first filtered based on their aggregated SMD
#' across the specified `smd_groups` and `smd_variables`, using the method
#' defined by `smd_type`. Then, among the remaining models, the best-performing
#' one(s) are selected based on the percentage of matched samples - either
#' overall or in the specified treatment groups (`perc_matched`).
#'
#' @return An S3 object of class `select_result`, containing the filtered and
#' prioritized optimization results. The object includes:
#'
#' - A `data.frame` with selected parameter combinations and performance
#' metrics.
#' - **Attribute `param_df`**: A `data.frame` with full parameter specifications
#'  (`iter_ID`, GPS/matching parameters, etc.), useful for manually refitting or
#'   reproducing results.
#'
#' The object also includes a custom `print()` method that summarizes:
#' - Number of selected combinations per SMD bin
#' - Corresponding aggregated SMD (mean or max)
#' - Overall or group-specific percentage matched
#'
#' @examples
#' \donttest{
#' # Define formula and set up optimization
#' formula_cancer <- formula(status ~ age * sex)
#' opt_args <- make_opt_args(cancer, formula_cancer, gps_method = "m1")
#' withr::with_seed(8252, {
#'   opt_results <- optimize_gps(
#'     data = cancer,
#'     formula = formula_cancer,
#'     opt_args = opt_args,
#'     n_iter = 2000
#'   )
#' })
#' # Select optimal combinations prioritizing SMD balance and matching in key
#' # groups
#' select_opt(
#'   x = opt_results,
#'   smd_groups = list(
#'     c("adenoma", "control"),
#'     c("control", "crc_beningn"),
#'     c("crc_malignant", "control")
#'   ),
#'   smd_variables = "age",
#'   smd_type = "max",
#'   perc_matched = c("adenoma", "crc_malignant")
#' )
#' }
#' @export
select_opt <- function(x,
                       smd_groups = NULL,
                       smd_variables = NULL,
                       smd_type = c("mean", "max"),
                       perc_matched = NULL) {
  # get valid treatments and their length
  treat_names <- attr(x, "treat_names")
  treat_length <- length(treat_names)

  # model_covs and their length
  model_covs <- attr(x, "model_covs")
  model_length <- length(model_covs)

  # validate smd_groups
  # The input should look like this:
  # smd_groups = list(
  #   c("GroupA", "GroupB"),
  #   c("GroupA", "GroupC")
  # )

  # for each pair
  if (!is.null(smd_groups)) {
    # ensure smd_groups is a list
    .chk_cond(
      !is.list(smd_groups) || !all(vapply(
        smd_groups,
        is.character,
        logical(1)
      )),
      "`smd_groups` must be a list of character vectors."
    )

    # for each pair
    for (i in seq_along(smd_groups)) {
      pair <- smd_groups[[i]]

      # only two items allowed
      .chk_cond(
        length(pair) != 2,
        sprintf("Each comparison in `smd_groups` must contain
              exactly two group names. Problem in element %d.", i)
      )

      # check for invalid values
      invalid <- pair[pair %nin% treat_names]
      .chk_cond(
        length(invalid) > 0,
        sprintf(
          "The following values in `smd_groups[[%d]]` are
        not valid group names: %s",
          i,
          paste(invalid, collapse = ", ")
        )
      )

      # check for duplicates
      .chk_cond(
        duplicated(pair)[1],
        sprintf("Comparison %d in `smd_groups` must contain
              two distinct group names.", i)
      )
    }

    # Normalize the smd_groups list --> sorting alphabetically
    smd_key <- lapply(smd_groups, function(g) sort(g))
    smd_key <- do.call(rbind, smd_key)
  } else {
    # smd_groups is NULL: create all unique pairwise combinations of treat_names
    combs <- t(utils::combn(unique(as.character(treat_names)), 2))

    # Normalize each pair by sorting alphabetically
    smd_key <- t(apply(combs, 1, function(x) sort(as.character(x))))
  }

  # data frame with all/selected pairwise comparisons for further reference
  smd_key_df <- data.frame(
    group1 = smd_key[, 1],
    group2 = smd_key[, 2],
    stringsAsFactors = FALSE
  )

  # validation function for smd_variables and perc_matched
  validate_arg <- function(x,
                           x_name,
                           available_values,
                           desired_length,
                           include_border = FALSE) {
    # Ensure it's a character vector (single or multiple)
    .chk_cond(
      !.check_vecl(x, NULL, FALSE),
      sprintf("The argument `%s` must be a single string
              or a vector of characters.", x_name)
    )

    chk::chk_character(x)

    # Check if all values are in the allowed set
    invalid <- x[x %nin% available_values]

    .chk_cond(
      length(invalid) > 0,
      sprintf(
        "The following values in `%s` are not allowed: %s.
        The set of allowed values is: %s",
        x_name,
        paste(invalid, collapse = ", "),
        word_list(available_values, quotes = TRUE)
      )
    )

    # Check for duplicates
    .chk_cond(
      any(duplicated(x)),
      sprintf("Duplicates are not allowed inside the `%s` argument.", x_name)
    )

    # Check length constraint
    if (include_border) {
      condition <- length(x) > desired_length
    } else {
      condition <- length(x) >= desired_length
    }

    .chk_cond(
      condition,
      sprintf(
        "The maximum allowed length of `%s` is %d. Setting `%s` to `NULL`.",
        x_name, desired_length - 1, x_name
      )
    )

    # Reset to NULL if too long
    if (length(x) >= desired_length) x <- NULL

    return(x)
  }

  ## process the smd_variables argument
  if (!is.null(smd_variables)) {
    smd_variables <- validate_arg(smd_variables,
      "smd_variables",
      model_covs,
      model_length,
      include_border = TRUE
    )
  }

  ## it has to be done after and separately, cause the upper check can reassign
  ## the value to NULL
  if (is.null(smd_variables)) {
    # use global smd max
    smd_variables <- model_covs
  }

  ## process the perc_matched argument
  if (!is.null(perc_matched)) {
    perc_matched <- validate_arg(
      perc_matched,
      "perc_matched",
      treat_names,
      treat_length
    )
  }

  ## the same thing as up
  if (is.null(perc_matched)) {
    # use global smd max
    perc_colnames <- "perc_matched"
  } else {
    # process to valid colname
    perc_colnames <- paste0("p_", perc_matched)
  }

  # Helper: normalize group pairs to alphabetical order
  normalize_pair <- function(g1, g2) {
    m <- mapply(function(a, b) sort(c(a, b)), g1, g2)
    data.frame(group1 = m[1, ], group2 = m[2, ], stringsAsFactors = FALSE)
  }

  ## process smd ===============================================================

  # extract and preprocess tha data
  smd_df <- attr(x, "smd_results")

  # take only complete cases (this mean cases where smd is defined)
  smd_df <- smd_df[stats::complete.cases(smd_df), ]

  # select the columns based on smd_variables
  smd_df <- smd_df[, c("iter_ID", "group1", "group2", smd_variables)]

  # Normalize the data frame group pairs
  normalized_df <- normalize_pair(smd_df$group1, smd_df$group2)
  smd_df$group1 <- normalized_df$group1
  smd_df$group2 <- normalized_df$group2

  # select the rows based on smd_groups
  # join with smd_key_df to filter out rows
  smd_df <- merge(smd_df,
    smd_key_df,
    by = c("group1", "group2")
  )

  # calculate mean/max by iterID for them

  # get the numeric columns after iter_ID
  metric_cols <- setdiff(names(smd_df), c("group1", "group2", "iter_ID"))

  # row-wise
  # Aggregate by group
  agg_fun <- if (smd_type == "mean") mean else max
  agg_df <- stats::aggregate(smd_df[metric_cols],
    by = smd_df["iter_ID"],
    FUN = agg_fun
  )

  # column-wise
  # Compute row-wise mean or max across metric columns
  metric_mat <- as.matrix(agg_df[metric_cols])
  row_stat <- if (smd_type == "mean") {
    rowMeans(metric_mat)
  } else {
    apply(metric_mat, 1, max)
  }

  # Append row-wise stat
  agg_df$overall_stat <- row_stat

  # categorize smd's using cut
  # Cut into SMD categories
  breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, Inf)
  labels <- c(
    "0-0.05", "0.05-0.10", "0.10-0.15", "0.15-0.20",
    "0.20-0.25", "0.25-0.30", "0.30-0.35", "0.35-0.40",
    "0.40-0.45", "0.45-0.50", ">0.50"
  )

  # assembling the output
  agg_df$smd_group <- cut(agg_df$overall_stat,
    breaks = breaks,
    labels = labels,
    right = FALSE
  )

  ## process percent_matched =================================================
  perc_df <- attr(x, "opt_results")

  ## select desired columns
  perc_df <- perc_df[, c("iter_ID", perc_colnames), drop = FALSE]

  ## filter out NAs
  perc_df <- perc_df[stats::complete.cases(perc_df), ]

  perc_components <- NULL # will store which columns were used

  ## if single column after iter_ID, then no need to aggregate
  if (ncol(perc_df) == 2L) {
    if (colnames(perc_df)[2L] != "perc_matched") {
      colnames(perc_df)[2L] <- "perc_matched"
    }
    perc_components <- colnames(perc_df)[2L]

    final_df <- merge(
      agg_df,
      perc_df[, c("iter_ID", "perc_matched"), drop = FALSE],
      by = "iter_ID"
    )
  } else {
    ## if more than two columns we take the *average* percentage
    metric_cols <- setdiff(names(perc_df), "iter_ID")
    perc_components <- metric_cols

    perc_df$perc_matched <- rowMeans(perc_df[metric_cols], na.rm = TRUE)

    final_df <- merge(
      agg_df,
      perc_df[, c("iter_ID", "perc_matched"), drop = FALSE],
      by = "iter_ID"
    )
  }

  # filter out NAs and Infs
  final_df <- final_df[apply(
    final_df, 1,
    function(row) all(!is.na(row) & !is.infinite(row))
  ), ]

  ## now filter to get the best rows
  # Get unique SMD groups
  groups <- unique(final_df$smd_group)

  # Initialize empty list to collect best rows
  best_rows_list <- list()

  # Loop through each smd group and extract rows with min perc_matched
  for (g in groups) {
    group_rows <- final_df[final_df$smd_group == g, ]

    # choose max perc_matched
    max_value <- max(group_rows$perc_matched, na.rm = TRUE)
    best_rows <- group_rows[group_rows$perc_matched == max_value, ]

    # choose min overall_stat
    min_value <- min(best_rows$overall_stat, na.rm = TRUE)
    best_rows <- best_rows[best_rows$overall_stat == min_value, ]

    best_rows_list[[as.character(g)]] <- best_rows
  }

  # Combine all best rows into one data.frame
  best_rows_final <- do.call(rbind, best_rows_list)

  ## reset row_names
  rownames(best_rows_final) <- NULL

  ## build param_df ============================================================
  full_df <- attr(x, "opt_results")

  remove_cols <- c(
    "perc_matched", "smd", "smd_group"
  )

  param_df <- merge(
    full_df[, colnames(full_df) %nin% remove_cols, drop = FALSE],
    best_rows_final[, colnames(best_rows_final) %nin% "perc_matched",
                    drop = FALSE],
    by = "iter_ID"
  )

  ## assemble structured result ===============================================

  # start with the data and desired classes
  res <- structure(
    best_rows_final,
    class = c("select_result", "data.frame")
  )

  # copy over relevant attributes from the original best_opt_result
  attrs_to_copy <- c(
    "treat_names",
    "model_covs",
    "optimization_time",
    "combinations_tested",
    "opt_results",
    "smd_results"
  )
  for (nm in attrs_to_copy) {
    attr(res, nm) <- attr(x, nm, exact = TRUE)
  }

  # add the parameter data.frame and function call
  attr(res, "param_df") <- param_df
  attr(res, "function_call") <- match.call()

  # printing is handled by print() methods
  return(res)
}

#' @export
summary.select_result <- function(object, digits = 3, ...) {
  x <- object

  # Remove rows with NA in key columns
  x <- x[!is.na(x$smd_group) &
    !is.na(x$overall_stat) &
    !is.na(x$perc_matched), , drop = FALSE]

  if (!nrow(x)) {
    # empty summary
    result_table <- data.frame(
      smd_group = character(0),
      unique_configs = integer(0),
      smd = numeric(0),
      perc_matched = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    # Split by group
    groups <- split(x, x$smd_group)

    # Keep only groups with valid overall_stat values
    groups <- groups[vapply(
      groups,
      function(g) any(!is.na(g$overall_stat)),
      logical(1L)
    )]

    # Sort groups by first valid overall_stat
    smd_order <- vapply(
      groups,
      function(g) min(g$overall_stat, na.rm = TRUE),
      numeric(1L)
    )
    sorted_groups <- names(sort(smd_order))

    # Build result table (unrounded here; rounding is in print())
    result_table <- data.frame(
      smd_group = character(),
      unique_configs = integer(),
      smd = numeric(),
      perc_matched = numeric(),
      stringsAsFactors = FALSE
    )

    for (grp in sorted_groups) {
      rows <- groups[[grp]]
      unique_configs <- length(unique(rows$iter_ID))
      first_row <- rows[1L, ]

      result_table <- rbind(
        result_table,
        data.frame(
          smd_group = grp,
          unique_configs = unique_configs,
          smd = first_row$overall_stat,
          perc_matched = first_row$perc_matched,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  res <- structure(
    result_table,
    digits          = as.integer(digits),
    perc_components = attr(object, "perc_components", exact = TRUE),
    class           = c("summary.select_result", "data.frame")
  )

  res
}

#' @export
print.summary.select_result <- function(x,
                                        digits = attr(x, "digits",
                                          exact = TRUE
                                        ) %||% 3L,
                                        ...) {
  digits <- as.integer(digits)
  if (!is.finite(digits) || digits < 0L) digits <- 3L

  cat("\n")
  cat("Optimization Selection Summary\n")
  cat("====================\n\n")

  # Info about what perc_matched represents
  comps <- attr(x, "perc_components", exact = TRUE)
  if (!is.null(comps)) {
    if (identical(comps, "perc_matched")) {
      cat("perc_matched: global percentage of matched units.\n\n")
    } else {
      cat(
        "perc_matched: mean percentage matched across groups: ",
        paste(comps, collapse = ", "),
        ".\n\n",
        sep = ""
      )
    }
  }

  result_table <- x

  # Table formatting
  header <- c("smd_group", "unique_configs", "smd", "perc_matched")
  col_widths <- c(12, 16, 8, 14)
  total_width <- sum(col_widths) + length(col_widths) + 1

  cat(strrep("-", total_width), "\n")
  cat(
    "|",
    paste(mapply(function(name, width) {
      format(name, width = width, justify = "centre")
    }, header, col_widths), collapse = "|"),
    "|\n",
    sep = ""
  )
  cat(strrep("-", total_width), "\n")

  if (nrow(result_table) > 0L) {
    for (i in seq_len(nrow(result_table))) {
      row <- result_table[i, ]

      cat(
        "|",
        format(row$smd_group, width = col_widths[1], justify = "left"), "|",
        format(row$unique_configs, width = col_widths[2], justify = "right"),
        "|",
        format(round(row$smd, digits),
          width = col_widths[3], justify = "right"
        ), "|",
        format(round(row$perc_matched, digits),
          width = col_widths[4], justify = "right"
        ),
        "|\n",
        sep = ""
      )
    }
  } else {
    cat(
      "|",
      format("<no rows>", width = total_width - 2L, justify = "centre"),
      "|\n",
      sep = ""
    )
  }

  cat(strrep("-", total_width), "\n\n")

  invisible(x)
}


#' @export
print.select_result <- function(x, n_show = 10L, ...) {
  cli::cli_text("{.strong select_result object} (selected matching
                configurations)")
  cli::cli_text("Subclass of {.cls best_opt_result} and {.cls data.frame}.")
  cli::cli_text("")

  # Optional info about what perc_matched represents
  comps <- attr(x, "perc_components", exact = TRUE)
  if (!is.null(comps)) {
    if (identical(comps, "perc_matched")) {
      cli::cli_text("perc_matched: global percentage of matched units.")
    } else {
      cli::cli_text(
        "perc_matched: mean percentage matched across groups:
        {.field {paste(comps, collapse = ', ')}}"
      )
    }
    cli::cli_text("")
  }

  # Reuse the core printer for best_opt_result if available
  if (exists(".print_best_opt_result_core", mode = "function")) {
    .print_best_opt_result_core(x, n_show = n_show, ...)
  } else {
    # Fallback: simple data.frame print if helper is not in scope
    n <- nrow(x)
    n_show <- as.integer(n_show)
    if (!is.finite(n_show) || n_show <= 0L) n_show <- 10L

    if (n > n_show) {
      cli::cli_text("Showing the first {n_show} of {n} rows:")
    } else {
      cli::cli_text("Showing all {n} rows:")
    }
    cli::cli_text("")
    base::print.data.frame(utils::head(x, n_show), ...)
  }

  invisible(x)
}

#' @export
str.select_result <- function(object, ...) {
  n <- nrow(object)
  p <- ncol(object)

  opt_time <- attr(object, "optimization_time", exact = TRUE)
  comb_test <- attr(object, "combinations_tested", exact = TRUE)
  opt_res <- attr(object, "opt_results", exact = TRUE)
  smd_res <- attr(object, "smd_results", exact = TRUE)
  treat_attr <- attr(object, "treat_names", exact = TRUE)
  model_covs <- attr(object, "model_covs", exact = TRUE)
  param_df <- attr(object, "param_df", exact = TRUE)
  comps <- attr(object, "perc_components", exact = TRUE)
  call <- attr(object, "function_call", exact = TRUE)

  ## Treatments -----------------------------------------------------------
  if (!is.null(treat_attr)) {
    if (is.factor(treat_attr)) {
      treat_levels <- levels(treat_attr)
    } else {
      treat_levels <- unique(as.character(treat_attr))
    }
  } else {
    treat_levels <- character(0)
  }
  k_treat <- length(treat_levels)

  ## SMD groups present in the selection --------------------------------
  smd_groups <- unique(object[["smd_group"]])
  n_smd_groups <- sum(!is.na(smd_groups))

  ## Header --------------------------------------------------------------
  cat("select_result object: selected GPS matching configurations\n")
  cat(sprintf(" Dimensions: %d rows x %d columns\n", n, p))

  cat(sprintf(
    " Optimization time (sec): %s\n",
    if (!is.null(opt_time)) as.character(opt_time) else "<unknown>"
  ))
  cat(sprintf(
    " Combinations tested in full optimization: %s\n",
    if (!is.null(comb_test)) as.character(comb_test) else "<unknown>"
  ))

  cat(sprintf(
    " Number of treatment levels: %d\n",
    k_treat
  ))
  cat(sprintf(
    " Treatment levels: %s\n",
    if (k_treat) paste(treat_levels, collapse = ", ") else "<none>"
  ))

  cat(sprintf(
    " Number of covariates in balance check: %d\n",
    length(model_covs %||% character(0))
  ))
  if (!is.null(model_covs)) {
    cat(sprintf(
      " Covariates: %s\n",
      paste(model_covs, collapse = ", ")
    ))
  }

  cat(sprintf(
    " SMD groups present in selection: %d\n",
    n_smd_groups
  ))
  if (n_smd_groups) {
    cat(sprintf(
      " SMD groups: %s\n",
      paste(stats::na.omit(as.character(smd_groups)), collapse = ", ")
    ))
  }

  if (!is.null(comps)) {
    if (identical(comps, "perc_matched")) {
      cat(" perc_matched: global percentage of matched units.\n")
    } else {
      cat(sprintf(
        " perc_matched: mean percentage matched across groups: %s\n",
        paste(comps, collapse = ", ")
      ))
    }
  }

  if (!is.null(opt_res)) {
    cat(sprintf(
      " opt_results (full grid): data.frame with %d rows and %d columns\n",
      NROW(opt_res), NCOL(opt_res)
    ))
  } else {
    cat(" opt_results: <none>\n")
  }

  if (!is.null(smd_res)) {
    cat(sprintf(
      " smd_results (full grid): data.frame with %d rows and %d columns\n",
      NROW(smd_res), NCOL(smd_res)
    ))
  } else {
    cat(" smd_results: <none>\n")
  }

  if (!is.null(param_df)) {
    cat(sprintf(
      " param_df (matched parameter rows): data.frame with %d rows and %d
      columns\n",
      NROW(param_df), NCOL(param_df)
    ))
  } else {
    cat(" param_df: <none>\n")
  }

  if (!is.null(call)) {
    cat(" Call:\n")
    cat(
      "  ",
      paste(deparse(call, width.cutoff = 80L), collapse = "\n  "),
      "\n",
      sep = ""
    )
  }

  cat("\nUnderlying data.frame structure:\n")

  ## Delegate to data.frame method for the actual structure -------------
  utils::str(unclass(object), ...)

  invisible(object)
}

#' @title Extract Parameter Grid for Selected Configurations
#'
#' @param x A \code{select_result} object returned by \code{select_opt()}.
#' @param smd_group Optional character vector of SMD groups
#'   (e.g. \code{"0.10-0.15"}). If provided, the returned data.frame
#'   is subsetted to these groups. If \code{NULL}, all rows are returned.
#'
#' @return A data.frame with the parameter combinations corresponding
#'   to the selected configurations.
#' @examples
#' \donttest{
#' # Define formula and set up optimization
#' formula_cancer <- formula(status ~ age * sex)
#' opt_args <- make_opt_args(cancer, formula_cancer, gps_method = "m1")
#' withr::with_seed(8252, {
#'   opt_results <- optimize_gps(
#'     data = cancer,
#'     formula = formula_cancer,
#'     opt_args = opt_args,
#'     n_iter = 2000
#'   )
#' })
#' # Select optimal combinations prioritizing SMD balance and matching in key
#' # groups
#' select_results <- select_opt(
#'   x = opt_results,
#'   smd_groups = list(
#'     c("adenoma", "control"),
#'     c("control", "crc_beningn"),
#'     c("crc_malignant", "control")
#'   ),
#'   smd_variables = "age",
#'   smd_type = "max",
#'   perc_matched = c("adenoma", "crc_malignant")
#' )
#'
#' # Extract the parameter grid from select_results for smd_group = "0.05-0.10"
#' get_select_params(select_results, smd_group = "0.05-0.10")
#' }
#' @export
get_select_params <- function(x, smd_group = NULL) {
  # Basic type check
  .chk_cond(
    !inherits(x, "select_result"),
    "`x` must be an object of class 'select_result'."
  )

  param_df <- attr(x, "param_df", exact = TRUE)

  .chk_cond(
    is.null(param_df),
    "Attribute 'param_df' is missing on this object."
  )

  if (!is.null(smd_group)) {
    # Ensure character vector
    smd_group <- as.character(smd_group)

    # Check that the column exists
    .chk_cond(
      !"smd_group" %in% colnames(param_df),
      "Column 'smd_group' is not present in 'param_df'."
    )

    # Subset to requested SMD groups
    param_df <- param_df[param_df$smd_group %in% smd_group, , drop = FALSE]
  }

  return(param_df)
}

#' @title Rerun GPS Estimation and Matching for a Selected Configuration
#'
#' @param x A \code{select_result} object returned by \code{select_opt()}.
#' @param data Data frame used in the original optimization (pass it the same
#'   way as in your original analysis, e.g. \code{data = cancer}).
#' @param formula Model formula used for GPS estimation (e.g.
#' \code{formula_cancer}).
#' @param smd_group Optional SMD bin to filter on
#'   (e.g. \code{"0.10-0.15"}). If \code{NULL}, no filtering is applied.
#' @param row Integer index of the row (after optional filtering by
#'   \code{smd_group}) to use. Defaults to \code{1}.
#' @param ... Extra args forwarded to \code{match_gps()}.
#'
#' @return The result of \code{match_gps()}.
#' @examples
#' \donttest{
#' # Define formula and set up optimization
#' formula_cancer <- formula(status ~ age * sex)
#' opt_args <- make_opt_args(cancer, formula_cancer, gps_method = "m1")
#' withr::with_seed(8252, {
#'   opt_results <- optimize_gps(
#'     data = cancer,
#'     formula = formula_cancer,
#'     opt_args = opt_args,
#'     n_iter = 2000
#'   )
#' })
#' # Select optimal combinations prioritizing SMD balance and matching in key
#' # groups
#' select_results <- select_opt(
#'   x = opt_results,
#'   smd_groups = list(
#'     c("adenoma", "control"),
#'     c("control", "crc_beningn"),
#'     c("crc_malignant", "control")
#'   ),
#'   smd_variables = "age",
#'   smd_type = "max",
#'   perc_matched = c("adenoma", "crc_malignant")
#' )
#'
#' # Extract the parameter grid from select_results for smd_group = "0.05-0.10"
#' get_select_params(select_results, smd_group = "0.05-0.10")
#'
#' # Rerun the analysis
#'
#' }
#' @export
run_selected_matching <- function(x,
                                  data,
                                  formula,
                                  smd_group = NULL,
                                  row = 1L,
                                  ...) {
  # sanity checks
  .chk_cond(
    !inherits(x, "select_result"),
    "`x` must be an object of class 'select_result'."
  )

  row <- as.integer(row)
  .chk_cond(
    !is.finite(row) || row < 1L,
    "`row` must be a positive integer."
  )

  # get parameter grid
  param_df <- get_select_params(x, smd_group = smd_group)

  .chk_cond(
    nrow(param_df) == 0L,
    "No parameter rows found for the requested `smd_group`."
  )

  .chk_cond(
    row > nrow(param_df),
    sprintf(
      "`row` (= %d) exceeds the number of available rows (= %d).",
      row, nrow(param_df)
    )
  )

  # take the chosen row
  par_row <- param_df[row, , drop = FALSE]

  ## --------- GPS estimation (via reconstructed call) ----------------------
  method_gps <- as.character(par_row[["method_gps"]])
  link <- as.character(par_row[["link"]])
  reference <- as.character(par_row[["reference"]])

  .chk_cond(
    any(is.na(c(method_gps, link, reference))),
    "Selected parameter row is missing `method_gps`, `link` or `reference`."
  )

  # use the original symbols from the user's call: data = cancer,
  # formula = formula_cancer, ...
  mc <- match.call(expand.dots = FALSE)

  est_call <- as.call(list(
    as.name("estimate_gps"),
    formula   = mc$formula,
    data      = mc$data,
    method    = method_gps,
    link      = link,
    reference = reference
  ))

  # evaluate in parent frame so that 'cancer', 'formula_cancer', etc. are found
  gps_mat <- eval(est_call, envir = parent.frame())

  # CSR with refitting
  csr_mat <- csregion(gps_mat)

  ## --------- Matching ----------------------------------------------------
  method_match <- as.character(par_row[["method_match"]])
  caliper <- as.numeric(par_row[["caliper"]])
  order <- as.character(par_row[["order"]])
  kmeans_cluster <- par_row[["kmeans_cluster"]]
  replace <- par_row[["replace"]]
  ties <- par_row[["ties"]]
  ratio <- par_row[["ratio"]]
  min_controls <- par_row[["min_controls"]]
  max_controls <- par_row[["max_controls"]]

  args_match <- list(
    csmatrix       = csr_mat,
    method         = method_match,
    caliper        = caliper,
    reference      = reference,
    order          = order,
    kmeans_cluster = kmeans_cluster,
    replace        = replace,
    ties           = ties,
    ratio          = ratio,
    min_controls   = min_controls,
    max_controls   = max_controls
  )

  # drop NA/NULL arguments
  args_match <- args_match[!vapply(
    args_match,
    function(z) {
      if (is.null(z)) {
        TRUE
      } else if (is.atomic(z) && length(z) == 1L) {
        is.na(z)
      } else {
        FALSE
      }
    },
    logical(1L)
  )]

  # add extra user args
  extra_args <- list(...)
  if (length(extra_args)) {
    args_match <- c(args_match, extra_args)
  }

  matched <- do.call(match_gps, args_match)

  # provenance
  attr(matched, "select_row") <- par_row
  attr(matched, "select_smd_group") <- par_row[["smd_group"]]
  attr(matched, "select_call") <- match.call()

  return(matched)
}

#' @title Define the Optimization Parameter Space for Matching
#'
#' @description `make_opt_args()` creates an object of class `"opt_args"` that
#' defines the parameter search space for `optimize_gps()`.
#'
#' The function accepts **vectors of values** for each customizable argument
#' involved in GPS estimation and matching. It computes the **Cartesian
#' product** of all parameter combinations, which serves as the input search
#' space for the random search algorithm used by `optimize_gps()`.
#'
#' To ensure valid optimization, the `data` and `formula` arguments must exactly
#' match those passed to `optimize_gps()`.
#'
#' @param data A `data.frame` containing all variables referenced in `formula`.
#'   Must match the dataset used in `optimize_gps()`.
#'
#' @param formula A valid formula specifying the treatment variable (left-hand
#'   side) and covariates (right-hand side). Interaction terms can be included
#'   using `*`. Must match the formula used in `optimize_gps()`.
#'
#' @param reference A single string or vector of treatment group levels to be
#'   used as the reference (baseline) group in both GPS estimation and matching.
#'
#' @param gps_method A string or vector of strings specifying GPS estimation
#'   methods. Allowed values are `"m1"` to `"m10"`. See *Details* below.
#'
#' @param matching_method A string or vector of strings specifying the matching
#'   method(s) to evaluate. Currently supported options are `"nnm"` and
#'   `"fullopt"`. See [match_gps()].
#'
#' @param caliper A numeric value or vector of values specifying caliper widths
#'   (i.e., maximum allowed GPS distance for matching). Same as in
#'   [match_gps()], but allows multiple values.
#'
#' @param order A string or vector of strings indicating the sorting order of
#'   logit-transformed GPS values before matching. Options are:
#' - `"desc"`: sort from highest to lowest (default),
#' - `"asc"`: sort from lowest to highest,
#' - `"original"`: keep original order,
#' - `"random"`: randomize order (use [set.seed()] for reproducibility).
#'
#' @param cluster An integer or vector of integers specifying the number of
#'   clusters for k-means clustering (if applicable).
#'
#' @param replace Logical value or vector of logicals indicating whether to
#'   allow matching with replacement. Same meaning as in [match_gps()], but
#'   supports multiple settings.
#'
#' @param ties Logical value or vector of logicals defining how ties should be
#'   handled during nearest-neighbor matching.
#'
#' @param ratio A numeric value or vector specifying the ratio of control to
#'   treated units for matching (used in `"nnm"`).
#'
#' @param min_controls A scalar or vector specifying the **minimum** number of
#'   controls to be matched to each treated unit (used in `"fullopt"`).
#'
#' @param max_controls A scalar or vector specifying the **maximum** number of
#'   controls to be matched to each treated unit (used in `"fullopt"`).
#'
#' @details The returned object is of class `"opt_args"` and is intended to be
#' passed directly to `optimize_gps()`. Internally, the function calculates the
#' full Cartesian product of all supplied parameter values and validates the
#' structure of each.
#'
#' The `gps_method` argument must contain one or more of the following codes:
#'
#' ```
#' | gps_method |      Method      |       Link Function         |
#' |------------|------------------|-----------------------------|
#' |    "m1"    |    multinom      |   generalized_logit         |
#' |    "m2"    |     polr         |   logistic                  |
#' |    "m3"    |     polr         |   probit                    |
#' |    "m4"    |     polr         |   loglog                    |
#' |    "m5"    |     polr         |   cloglog                   |
#' |    "m6"    |     polr         |   cauchit                   |
#' |    "m7"    |     vglm         |   multinomial_logit         |
#' |    "m8"    |     vglm         |   reduced_rank_ml           |
#' |    "m9"    |    brglm2        |   baseline_category_logit   |
#' |   "m10"    |    mblogit       |   baseline_category_logit   |
#' ```
#'
#' The object includes a custom S3 `print()` method that displays:
#' - A summary table of all allowed values for each optimization parameter,
#' - The total number of unique parameter combinations (i.e., the size of the
#' search space).
#'
#' @return An object of class `"opt_args"`, containing all valid parameter
#'   combinations to be sampled by `optimize_gps()`. Use `print()` to explore
#'   the defined search space.
#'
#' @seealso [optimize_gps()], [match_gps()], [estimate_gps()]
#'
#' @examples
#' # Define formula and dataset
#' formula_cancer <- formula(status ~ age * sex)
#'
#' # Create search space with multiple values for GPS and matching
#' opt_args <- make_opt_args(
#'   data = cancer,
#'   formula = formula_cancer,
#'   gps_method = c("m1", "m2", "m9"),
#'   matching_method = c("nnm", "fullopt"),
#'   caliper = c(0.1, 0.2),
#'   order = c("desc", "random"),
#'   reference = "control"
#' )
#'
#' # Print summary of the search space
#' print(opt_args)
#'
#' @export
make_opt_args <- function(
  data = NULL,
  formula,
  reference = NULL,
  gps_method = paste0("m", 1:10),
  matching_method = c("fullopt", "nnm"),
  caliper = seq(0.01, 10, 0.01),
  order = c("desc", "asc", "original", "random"),
  cluster = 2,
  replace = c(TRUE, FALSE),
  ties = c(TRUE, FALSE),
  ratio = 1,
  min_controls = 1,
  max_controls = 1
) {
  ############## PARAMETER CHECK ###############################################
  ## gps_method
  allowed_gps_methods <- paste0("m", 1:10)
  validate_optarg(
    gps_method,
    allowed_gps_methods
  )

  ## reference
  if (!is.null(data)) .check_df(data)

  # formula
  data_list <- .process_formula(formula, data)

  # reference processing
  treatment_var <- data_list[["treat"]]
  treatment_levels <- as.character(unique(treatment_var))

  # reference
  if (is.null(reference)) reference <- treatment_levels[1]

  ref_list <- lapply(reference, function(ref) {
    .process_ref(
      treatment_var,
      ordinal_treat = NULL,
      reference = ref
    )[["reference"]]
  })
  reference <- unlist(ref_list)

  validate_optarg(
    reference,
    treatment_levels
  )

  ## matching_method
  allowed_matching_methods <- c("fullopt", "nnm")
  validate_optarg(
    matching_method,
    allowed_matching_methods
  )

  ## caliper
  .chk_cond(
    !.check_vecl(caliper),
    "The argument `caliper` has to be a numeric vector."
  )
  .chk_cond(
    any(duplicated(caliper)),
    "Duplicates inside the `caliper` vector are not allowed."
  )
  .chk_cond(
    any(caliper <= 0),
    "Zeros and negative values are not allowed inside the `caliper` vector."
  )

  ## order
  allowed_orders <- c("desc", "asc", "original", "random")
  validate_optarg(
    order,
    allowed_orders
  )

  ## cluster
  allowed_clusters <- seq_along(treatment_levels)
  validate_optarg(
    cluster,
    allowed_clusters,
    quotes = FALSE
  )

  ## replace, ties and ratio
  if ("nnm" %in% matching_method) {
    allowed_ratio <- 1:20
    validate_optarg(
      ratio,
      allowed_ratio,
      quotes = FALSE
    )

    allowed_replace <- c(TRUE, FALSE)
    validate_optarg(
      replace,
      allowed_replace
    )

    allowed_ties <- c(TRUE, FALSE)
    validate_optarg(
      ties,
      allowed_ties
    )
  }

  ## max and min_controls
  if ("fullopt" %in% matching_method) {
    allowed_controls <- 1:20

    validate_optarg(
      min_controls,
      allowed_controls,
      quotes = FALSE
    )

    validate_optarg(
      max_controls,
      allowed_controls,
      quotes = FALSE
    )

    ## here you could still enforce max_controls >= min_controls if needed
  }

  # remove variable arguments if only one method specified
  if (length(matching_method) == 1 && matching_method == "nnm") {
    min_controls <- NULL
    max_controls <- NULL
  } else if (length(matching_method) == 1 && matching_method == "fullopt") {
    replace <- NULL
    ties <- NULL
    ratio <- NULL
  }

  ############# DEFINE THE LIST ################################################
  user_args <- list(
    gps_method      = gps_method,
    reference       = reference,
    matching_method = matching_method,
    caliper         = caliper,
    order           = order,
    cluster         = cluster,
    replace         = replace,
    ties            = ties,
    ratio           = ratio,
    min_controls    = min_controls,
    max_controls    = max_controls
  )

  # drop NULL components
  user_args <- user_args[!vapply(user_args, is.null, logical(1L))]

  ############# CALCULATE NUMBER OF COMBINATIONS ###############################
  if (length(unique(user_args[["matching_method"]])) == 2L) {
    total_combinations <- 0L

    for (method in user_args[["matching_method"]]) {
      method_args <- user_args

      if (method == "nnm") {
        method_args[c("min_controls", "max_controls")] <- NULL
      } else if (method == "fullopt") {
        method_args[c("replace", "ties", "ratio")] <- NULL
      } else {
        warning(sprintf("Unknown matching method: '%s'. Skipping.", method))
        next
      }

      method_args$matching_method <- method

      unique_lengths <- vapply(
        method_args,
        function(x) length(unique(x)),
        integer(1L)
      )
      total_combinations <- total_combinations + prod(unique_lengths)
    }
  } else {
    unique_lengths <- vapply(
      user_args,
      function(x) length(unique(x)),
      integer(1L)
    )
    total_combinations <- prod(unique_lengths)
  }

  total_combinations_fmt <- format(total_combinations, scientific = FALSE)
  model_covs <- colnames(data_list[["model_covs"]])

  ############# RETURN STRUCTURED OBJECT ######################################
  res <- structure(
    user_args,
    total_combinations = total_combinations_fmt,
    model_covs         = model_covs,
    class              = c("opt_args", "list")
  )

  return(res)
}


#' @export
## Custom print method for the S3 class
print.opt_args <- function(x, ...) {
  cat("Optimization Argument Set (class: opt_args)\n")
  cat(strrep("-", 40), "\n")

  # Print each argument and its unique values
  for (name in names(x)) {
    values <- unique(x[[name]])
    val_str <- if (length(values) > 10) {
      paste0("[", length(values), " values]")
    } else {
      paste(values, collapse = ", ")
    }
    cat(sprintf("%-16s: %s\n", name, val_str))
  }

  # Print total combinations
  total_combinations <- attr(x, "total_combinations")
  if (!is.null(total_combinations)) {
    cat(strrep("-", 40), "\n")
    cat("Total combinations:", total_combinations, "\n")
  }

  invisible(x)
}

#' @export
str.opt_args <- function(object, ...) {
  # pull attributes
  total_combinations <- attr(object, "total_combinations")
  model_covs <- attr(object, "model_covs")

  n_args <- length(object)
  n_covs <- if (!is.null(model_covs)) length(model_covs) else NA_integer_

  cat("Optimization argument grid (class: opt_args)\n")

  if (!is.null(total_combinations)) {
    cat(sprintf("  Total combinations in grid: %s\n", total_combinations))
  }

  cat(sprintf("  Number of argument dimensions: %d\n", n_args))

  if (!is.null(model_covs)) {
    cat(sprintf("  Number of model covariates: %d\n", n_covs))
    # show only first few covariates if there are many
    if (n_covs <= 10L) {
      cat("  Model covariates: ",
        paste(model_covs, collapse = ", "),
        "\n",
        sep = ""
      )
    } else {
      cat("  Model covariates: ",
        paste(utils::head(model_covs, 10L), collapse = ", "),
        ", ...\n",
        sep = ""
      )
    }
  }

  cat("\nArgument dimensions (unique values per argument):\n")

  # small per-argument summary (no full dump, str() will follow)
  for (name in names(object)) {
    values <- object[[name]]
    uvals <- unique(values)
    n_u <- length(uvals)
    cls <- paste(class(values), collapse = "/")

    # short preview of values
    preview <- if (n_u == 0L) {
      "<empty>"
    } else if (is.numeric(uvals)) {
      sprintf(
        "numeric, range [%.3g, %.3g], %d unique",
        min(uvals), max(uvals), n_u
      )
    } else if (is.logical(uvals)) {
      paste0("logical: ", paste(uvals, collapse = ", "), " (", n_u, " unique)")
    } else {
      if (n_u <= 5L) {
        paste0(
          cls, ": ",
          paste(uvals, collapse = ", "),
          " (", n_u, " unique)"
        )
      } else {
        paste0(
          cls, ": ",
          paste(utils::head(uvals, 5L), collapse = ", "),
          ", ... (", n_u, " unique)"
        )
      }
    }

    cat(sprintf("  $ %-15s: %s\n", name, preview))
  }

  cat("\nUnderlying list structure:\n")
  utils::str(unclass(object), ...)

  invisible(object)
}
