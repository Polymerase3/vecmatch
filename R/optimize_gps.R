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
#' @param n_cores Integer. Number of CPU cores to use for parallel execution. If
#'   set to a value greater than 1, a parallel backend is registered using
#'   [future::multisession()]. **Note:** parallel execution can significantly
#'   increase memory usage. With large datasets or high `n_iter` values, RAM
#'   consumption may spike, especially on systems with 16-32 GB RAM. Users are
#'   advised to monitor system resources. Internally, the function performs
#'   memory cleanup post-execution to manage resource usage efficiently.
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
#' - `method_match`: Matching method used in [match_gps()], e.g., `"nnm"` or `"fullopt"`
#' - `caliper`: Caliper value used in [match_gps()]
#' - `order`: Ordering of GPS scores prior to matching
#' - `kmeans_cluster`: Number of k-means clusters used
#' - `replace`: Whether replacement was used in matching (`nnm` only)
#' - `ties`: Tie-breaking rule in nearest-neighbor matching (`nnm` only)
#' - `ratio`: Control-to-treated ratio for `nnm`
#' - `min_controls`, `max_controls`: Minimum and maximum controls for `fullopt`
#' - `reference`: Reference group used in both [estimate_gps()] and [match_gps()]
#' - `perc_matched`: Percentage of matched samples (from [balqual()])
#' - `smd`: Maximum standardized mean difference (from [balqual()])
#' - `p_{group_name}`: Percent matched per treatment group (based on group sample size)
#' - `method_gps`: GPS estimation method used (from [estimate_gps()])
#' - `link`: Link function used in GPS model
#' - `smd_group`: SMD range category for the row
#'
#' The resulting `best_opt_result` object also includes a custom `print()`
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
#' - **`opt_results`**: A `data.frame` of optimization results. Each row corresponds to a unique parameter combination tested. For a complete description of columns, see the ***Details*** section.
#'
#' - **`optimization_time`**: Time (in seconds) taken by the optimization loop (i.e., the core `for`-loop that evaluates combinations). This does **not** include the time needed for GPS estimation, pre-processing, or merging of results after loop completion. On large datasets, these excluded steps can still be substantial.
#'
#' - **`combinations_tested`**: Total number of unique parameter combinations evaluated during optimization.
#'
#' - **`smd_results`**: A detailed table of standardized mean differences (SMDs) for all pairwise treatment group comparisons and for all covariates specified in the `formula`. This is used by the [select_opt()] function to filter optimal models based on covariate-level balance across groups.
#'
#' - **`treat_names`**: A character vector with the names of the unique treatment groups.
#'
#' - **`model_covs`**: A character vector listing the model covariates (main effects and interactions) used in the `formula`. These names correspond to the variables shown in the `smd_results` table.
#'
#' @examples
#' # Define formula for GPS estimation and matching
#' formula_cancer <- formula(status ~ age * sex)
#'
#' # Set up the optimization parameter space
#' opt_args <- make_opt_args(cancer, formula_cancer, gps_method = "m1")
#'
#' # Run optimization with 2000 random parameter sets and a fixed seed
#' \dontrun{
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
                         n_cores = 1,
                         opt_args = NULL) {
  ####################### INPUT CHECKING #######################################
  ########################### AND ##############################################
  ####################### DATA PROCESSING ######################################
  # defining a placeholder func for callr (function() inside callr::r throws
  # an error)
  placeholder_func <- function() {

  }

  # isolating the execution of the optimize_gps inside a separate R session -->
  # after closing it, the memory is freed completely. Small RSS may remain but
  # it frees after closing the parent process or restarting R
  callr::r(placeholder_func,
    {
      # check that the n_iter is integer and greater than 1
      .check_integer(n_iter)
      .chk_cond(
        n_iter < 1,
        "The number of iterations (`n_iter`) has to be at least 1."
      )

      # check that n_cores is numeric and at least 1
      .check_integer(n_cores)
      .chk_cond(
        n_cores < 1,
        "You must assign at least one core to the process."
      )

      # -------- parallel backend ---------------------------------------------
      ## process the number of cores (register backend)
      if (n_cores > 1) {
        # ensure required packages for parallel backend are installed
        rlang::check_installed(
          c(
            "future", "foreach", "doFuture",
            "parallel", "progressr", "doRNG"
          ),
          reason = "to register parallel backend"
        )

        # explicitly import dorng
        #' @importFrom doRNG %dorng%

        # detect the number of available CPU cores
        all_cores <- parallel::detectCores()

        # if requested cores exceed available cores, reduce n_cores
        if (n_cores > all_cores) {
          new_cores <- max(all_cores - 2, 1)
          warning(sprintf(
            "Requested number of cores (%d) exceeds available cores (%d).
        Setting n_cores to %d.",
            n_cores, all_cores, new_cores
          ))
          n_cores <- new_cores
        }

        ## define startup and close function for the parallel backend
        startup_par <- function(n_cores) {
          # set up parallel backend using `doFuture`
          future::plan(future::multisession, workers = n_cores, gc = TRUE)
          # works for every platform

          # register foreach to use `future`
          doFuture::registerDoFuture()
        }

        ## close function --> closes parallel backend and frees RAM
        shutdown_par <- function() {
          ## has to free the RAM completely --> no remaining allocs after run
          ## 1 – switch to sequential so new futures stay local
          future::plan(future::sequential)

          ## 2 – collect zombie forks on Unix
          if (.Platform$OS.type == "unix") {
            parallel::mccollect(wait = FALSE)
          }

          ## 3 – two explicit garbage collections:
          invisible(gc(full = TRUE, reset = TRUE))
          invisible(gc(full = TRUE))
        }
      } else {
        # only the foreach and progressr packages necessary
        rlang::check_installed(c("foreach", "progressr"),
          reason = "to execute for loops"
        )

        #' @importFrom foreach %do%
      }

      # define the infix operator for the forloops
      `%doparallel%` <- if (n_cores == 1) {
        `%do%`
      } else {
        `%dorng%`
      }

      ## adding default optimization parameters if not specified
      if (is.null(opt_args)) {
        opt_args <- make_opt_args(
          data = data,
          formula = formula
        )
      }

      # validate class of opt_args if provided
      .chk_cond(
        !inherits(opt_args, "opt_args"),
        "`opt_args` must be of class 'opt_args'.
            Use `make_opt_args()` to create it."
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
      estimate_methods$gps_method <- paste0("m", 1:10)
      ## --> see the documentation of make_opt_args()!

      # subset for gps_method specified in the opt_args
      estimate_methods <- estimate_methods[
        estimate_methods$gps_method %in% opt_args[["gps_method"]],
      ]

      # Create reference levels from unique treatment values
      available_refs <- data.frame(refs = opt_args[["reference"]])

      # Cartesian product of methods x reference levels
      estimate_space <- merge(estimate_methods,
        available_refs,
        by = NULL
      )

      # Adding unique names to each column
      estimate_space$row_name <- paste(rep("estimate", nrow(estimate_space)),
        seq_len(nrow(estimate_space)),
        sep = "_"
      )

      # remove "polr" from the extimate space if ordinal treat is null
      if (is.null(ordinal_treat)) {
        estimate_space <- estimate_space[estimate_space$method != "polr", ]
      }

      ##############################################################################
      ###### MATCH GPS

      ## defining number of iterations and searches
      n_iter_final <- n_iter # final number of iterations used for looping
      n_iter <- n_iter * 10 # number of iterations using for grid search
      # many are duplicates --> remove --> less samples

      # generate base parameter grid using random sampling
      withr::with_preserve_seed({
        search_matching <- data.frame(
          gps_model = sample(unique(estimate_space$row_name),
            n_iter,
            replace = TRUE
          ),
          method = sample(opt_args[["matching_method"]],
            n_iter,
            replace = TRUE
          ),
          caliper = sample(opt_args[["caliper"]],
            n_iter,
            replace = TRUE
          ),
          order = sample(opt_args[["order"]],
            n_iter,
            replace = TRUE
          ),
          kmeans_cluster = sample(opt_args[["cluster"]],
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
        search_matching$replace[is_nnm] <- sample(opt_args[["replace"]],
          n_nnm,
          replace = TRUE
        )
      })

      withr::with_preserve_seed({
        search_matching$ties[is_nnm] <- sample(opt_args[["ties"]],
          n_nnm,
          replace = TRUE
        )
      })

      withr::with_preserve_seed({
        search_matching$ratio[is_nnm] <- sample(opt_args[["ratio"]],
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

      # note: all sample() calls are wraped inside with_preserve_seed --> without
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
      names(search_matching)[which(names(search_matching) == "refs")] <- "reference"
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
            search_matching[sample(nrow(search_matching),
              n_iter_final,
              replace = FALSE
            ), ]
          })
        },
        error = function(e) {
          warning("Sampling the final number of rows without replacement failed.
      Consider expanding the parameter range.
      Falling back to sampling with replacement.")
          search_matching[sample(nrow(search_matching),
            n_iter_final,
            replace = TRUE
          ), ]
        }
      )

      ###################### FITTING ESTIMATE GPS #################################
      # logging
      rlang::inform("Initiating estimation of the GPS...\n")

      # starting parallel backend if n_cores > 1
      if (n_cores > 1 && n_cores <= 5) {
        rlang::inform("Registering parallel backend...\n")
        startup_par(n_cores)
      } else if (n_cores > 5) {
        # adding more than 5 cores would become inefficient, especially with larger
        # datasets
        n_cores_aux <- 5
        rlang::inform("Registering parallel backend...\n")
        startup_par(n_cores_aux)
      }

      # Enable global handler for the progress (once is enough)
      progressr::handlers("txtprogressbar")

      # to avoid seed leaks
      withr::with_preserve_seed({
        # printing out the progress bar
        progressr::with_progress({
          ## defining the loop length
          loop_seq <- seq_len(nrow(estimate_space))
          p <- progressr::progressor(along = loop_seq)

          ## looping through all estimate space combinations
          suppressMessages({
            estimate_results <- foreach::foreach(
              i = loop_seq,
              .packages = c("vecmatch", "withr"), # add any other needed packages
              # i needed to export all used objects to the workers in the parallel
              # setup
              .export = c(
                "data",
                "formula",
                "estimate_space",
                "ordinal_treat",
                "estimate_gps",
                "csregion",
                "p"
              ),
              .errorhandling = "pass"
            ) %doparallel% {
              # i know it may seems unnecessary to wrap it all inside a function and
              # call it at the end, but it has a purpose. During the debugging i
              # noticed that the workers somehow share a common environment and
              # managed to change the values of some variables during the run. It
              # was a total mess.
              # so the solution was to wrap the whole loop inside a function and
              # call it at the end. Firstly, it creates an isolated environment for
              # each iteration (function env), and the whole iteration is executed
              # only within this environment, which is then deleted

              run_iteration <- function(i) {
                # no logging from the loop
                withr::local_options(list(foreach.quiet = TRUE))

                # defining the current argument list
                arglist_loop <- list(
                  data = data,
                  formula = formula,
                  method = estimate_space[i, "method"],
                  link = estimate_space[i, "link"],
                  reference = as.character(estimate_space[i, "refs"]),
                  ordinal_treat = ordinal_treat
                )

                # estimating the gps
                # double caution - better safe than sorry
                withr::with_preserve_seed({
                  tryCatch(
                    {
                      loop_estimate <- do.call(estimate_gps, arglist_loop)
                    },
                    error = function(e) {
                      list(
                        error = TRUE,
                        message = sprintf(
                          "Error with method `%s`: %s",
                          arglist_loop$method,
                          conditionMessage(e)
                        )
                      )
                    }
                  )

                  # calculating csregion borders
                  invisible(
                    utils::capture.output(loop_estimate <- csregion(loop_estimate))
                  )
                })

                # Update progress
                p(sprintf("Running %d/%d", i, max(loop_seq)))

                # remove objects
                rm(list = setdiff(ls(), "loop_estimate"))

                # defining the output
                loop_estimate <- list(loop_estimate)
                names(loop_estimate) <- estimate_space[i, "row_name"]
                return(loop_estimate)
              }

              # call the function within isolated env (function env)
              run_iteration(i)
            }
          })
        })
      })

      ## unnesting the estimate_results list
      list_names <- vapply(estimate_results, function(x) names(x[1]), character(1))
      estimate_results <- lapply(estimate_results, function(x) x[[1]])

      # assigning unique row_names for the estimates
      names(estimate_results) <- list_names

      # shutting down parallel backend if n_cores > 1
      if (n_cores > 1) {
        rlang::inform("Freeing RAM...\n")
        shutdown_par()
      }

      ###################### MATCHING OPTIMIZER ####################################
      # The search grid is already genereated
      # logging
      rlang::inform("Initiating matching optimization...\n")

      # starting parallel backend if n_cores > 1
      if (n_cores > 1) {
        rlang::inform("Registering parallel backend...\n")
        startup_par(n_cores)
      }

      # logging
      rlang::inform("Optimizing matching...\n")

      ## starting time
      time_start <- Sys.time()

      ## preserving the seed (avoiding seed leaks)
      withr::with_preserve_seed({
        suppressMessages({
          progressr::with_progress({
            ## loop length and progress bar
            loop_seq <- seq_len(n_iter_final)
            throttle <- 100
            p <- progressr::progressor(along = seq_len(n_iter_final / throttle))

            # Precompute unique treatment names outside the foreach loop
            all_treatments <- unique(estimate_results[[1]]$treatment)
            treatment_cols <- paste0("p_", all_treatments)

            # Get the structure of smd_df beforehand
            smd_colnames <- c(c("group1", "group2"), attr(opt_args, "model_covs"))

            # Predefine fallback
            smd_template <- as.data.frame(
              matrix(NA,
                nrow = 1,
                ncol = length(smd_colnames)
              ),
              stringsAsFactors = FALSE
            )
            colnames(smd_template) <- smd_colnames

            opt_results <- foreach::foreach(
              i = loop_seq,
              # i know it may seem weird but it actually works. It doesn't work with
              # the quick load from devtools, the package must be installed!
              # then all functions from the vecmatch are exported to the workers
              .packages = "vecmatch",

              # exporting all necessary objects
              .export = c(
                "search_matching", "estimate_results", "formula",
                "treatment_cols", "smd_colnames", "smd_template",
                "match_gps", "balqual"
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
                args_loop <- args_loop[-1]

                # Defaults in case of error
                perc_matched <- NA_real_
                smd <- NA_real_
                perc <- as.data.frame(t(rep(NA_real_, length(treatment_cols))))
                colnames(perc) <- treatment_cols

                smd_df <- smd_template # predefine smd_df with correct structure

                try(
                  {
                    # matching
                    loop_estimate <- do.call(match_gps, args_loop)

                    # max SMD and %matched statistics
                    utils::capture.output({
                      qual_out <- balqual(loop_estimate,
                        formula,
                        type = "smd",
                        statistic = "max",
                        round = 8,
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
                    computed <- stats::setNames(as.list(ptab$p), paste0(
                      "p_",
                      ptab$Treatment
                    ))

                    for (col in names(computed)) {
                      if (col %in% treatment_cols) {
                        perc[[col]] <- computed[[col]]
                      }
                    }

                    # update template with correct number of rows
                    smd_df <- as.data.frame(
                      matrix(NA,
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

                # Update progress
                if (i %% throttle == 0) {
                  p(sprintf(
                    "Running %d/%d",
                    i,
                    max(loop_seq)
                  ))
                }

                # setting up the resultin data frame
                result_row <- cbind(
                  iter_ID = iter_ID,
                  search_matching[i, ],
                  perc_matched = perc_matched,
                  smd = smd,
                  perc
                )

                # corresponding smd data.frame
                smd_df <- cbind(iter_ID = iter_ID, smd_df)

                # returning output from single iteration
                res <- list(
                  results = as.data.frame(result_row, stringsAsFactors = FALSE),
                  smd_dfs = as.data.frame(smd_df, stringsAsFactors = FALSE)
                )

                ## --- tidy-up: remove *everything* except 'result' ----------------
                rm(list = setdiff(ls(), "res"))

                ## optional: a *light* sweep every 256th iteration
                # if (i %% 1024 == 0L) gc(FALSE)

                # returnig the results
                return(res)
              }

              # running the function inside an isolated env
              run_iteration(i)
            }
          })
        })
      })

      # defining the stop time and running time
      time_stop <- Sys.time()
      time_diff <- round(as.numeric(time_stop - time_start, units = "secs"), 2)

      # shutting down parallel backend if n_cores > 1
      if (n_cores > 1) {
        rlang::inform("Freeing RAM...\n")
        shutdown_par()
      }

      ## extract looping results
      opt_results_df <- do.call(rbind, lapply(opt_results, `[[`, "results"))
      smd_df_all <- do.call(rbind, lapply(opt_results, `[[`, "smd_dfs"))

      # processing the results
      ## merge the gps_model
      opt_results_df <- merge(opt_results_df,
        estimate_space[, c("method", "link", "row_name")],
        by.x = "gps_model",
        by.y = "row_name"
      )

      ## rename method.y anx x
      names(opt_results_df)[
        which(names(opt_results_df) == "method.y")
      ] <- "method_gps"
      names(opt_results_df)[
        which(names(opt_results_df) == "method.x")
      ] <- "method_match"

      # Define SMD intervals
      breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, Inf)
      labels <- c(
        "0-0.05", "0.05-0.10", "0.10-0.15", "0.15-0.20",
        "0.20-0.25", "0.25-0.30", ">0.30"
      )

      # Cut SMD into intervals
      opt_results_df$smd_group <- cut(opt_results_df$smd,
        breaks = breaks,
        labels = labels,
        right = FALSE
      )

      ### BEST RESULTS OUTPUT FILTERING SECTION ####################################

      # Remove rows with missing perc_matched or smd_group
      filtered <- opt_results_df[
        !is.na(opt_results_df$perc_matched) & !is.na(opt_results_df$smd_group),
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

      # define custom s3 class
      class(best_rows_final) <- c("best_opt_result", class(best_rows_final))

      # add time_elapsed and total number of combinations tested
      attr(best_rows_final, "optimization_time") <- time_diff
      attr(best_rows_final, "combinations_tested") <- format(n_iter_final,
        scientific = FALSE
      )
      attr(best_rows_final, "opt_results") <- opt_results_df
      attr(best_rows_final, "smd_results") <- smd_df_all
      attr(best_rows_final, "treat_names") <- unique(
        estimate_results[[1]]$treatment
      )
      attr(best_rows_final, "model_covs") <- attr(opt_args, "model_covs")

      # printing out to the console
      rlang::inform("=========================================================\n\n")

      print(best_rows_final)
      return(invisible(best_rows_final))
    },
    spinner = TRUE
  )
}

#' @export
print.best_opt_result <- function(x, digits = 3, ...) {
  cat("Best Optimization Results by SMD Group\n")
  cat("======================================\n\n")

  # Filter and sort
  x <- x[!is.na(x$smd) & !is.na(x$perc_matched) & !is.na(x$smd_group), ]
  groups <- split(x, x$smd_group)
  smd_order <- sapply(groups, function(g) min(g$smd, na.rm = TRUE))
  sorted_groups <- names(sort(smd_order))

  # Build summary
  result_table <- data.frame(
    smd_group = character(),
    unique_configs = integer(),
    smd = numeric(),
    perc_matched = numeric(),
    stringsAsFactors = FALSE
  )

  for (grp in sorted_groups) {
    rows <- groups[[grp]]
    config_cols <- setdiff(names(rows), c("smd", "perc_matched"))
    unique_configs <- nrow(unique(rows[, config_cols]))
    smd_vals <- round(rows$smd, digits)
    perc_vals <- round(rows$perc_matched, digits)

    summary_rows <- unique(data.frame(
      smd_group = grp,
      unique_configs = unique_configs,
      smd = smd_vals,
      perc_matched = perc_vals
    ))

    result_table <- rbind(result_table, summary_rows)
  }

  # Table formatting
  header <- c("smd_group", "unique_configs", "smd", "perc_matched")
  col_widths <- c(12, 16, 8, 14)
  total_width <- sum(col_widths) + length(col_widths) + 1

  cat(strrep("-", total_width), "\n")
  cat("|", paste(mapply(function(name, width) {
    format(name, width = width, justify = "centre")
  }, header, col_widths), collapse = "|"), "|\n", sep = "")
  cat(strrep("-", total_width), "\n")

  for (i in seq_len(nrow(result_table))) {
    row <- result_table[i, ]
    cat("|",
      format(row$smd_group, width = col_widths[1], justify = "left"), "|",
      format(row$unique_configs,
        width = col_widths[2],
        justify = "right"
      ), "|",
      format(round(row$smd, digits),
        width = col_widths[3],
        justify = "right"
      ), "|",
      format(round(row$perc_matched, digits),
        width = col_widths[4],
        justify = "right"
      ), "|\n",
      sep = ""
    )
  }
  cat(strrep("-", total_width), "\n\n")

  # Optimization summary
  time_taken <- attr(x, "optimization_time")
  n_combos <- attr(x, "combinations_tested")

  if (!is.null(time_taken) && !is.null(n_combos)) {
    cat("Optimization Summary\n")
    cat("--------------------\n")
    cat("Total combinations tested  : ", format(n_combos, big.mark = ""),
      "\n",
      sep = ""
    )
    cat("Total optimization time [s]: ", format(time_taken), "\n", sep = "")
  }

  invisible(x)
}
#' @title Select Optimal Parameter Combinations from Optimization Results
#'
#' @description `select_opt()` is a helper function to filter and prioritize
#' results from `optimize_gps()` based on the specific goals of a study.
#' Depending on the research design, certain pairwise comparisons or treatment
#' groups may be more important than others. For example:
#'
#' - You may want to prioritize matching between a specific groups (e.g. specific disease vs. controls), while ignoring other group comparisons during SMD evaluation.
#' - You may wish to retain as many samples as possible from a critical group or set of groups, regardless of matching rates in other groups.
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
#' - A `data.frame` with selected parameter combinations and performance metrics.
#' - **Attribute `param_df`**: A `data.frame` with full parameter specifications (`iter_ID`, GPS/matching parameters, etc.), useful for manually refitting or reproducing results.
#'
#' The object also includes a custom `print()` method that summarizes:
#' - Number of selected combinations per SMD bin
#' - Corresponding aggregated SMD (mean or max)
#' - Overall or group-specific percentage matched
#'
#' @examples
#' # Define formula and set up optimization
#' formula_cancer <- formula(status ~ age * sex)
#' opt_args <- make_opt_args(cancer, formula_cancer, gps_method = "m1")
#' \dontrun{
#' withr::with_seed(8252, {
#'   opt_results <- optimize_gps(
#'     data = cancer,
#'     formula = formula_cancer,
#'     opt_args = opt_args,
#'     n_iter = 2000
#'   )
#' })
#' }
#' # Select optimal combinations prioritizing SMD balance and matching in key
#' # groups
#' \dontrun{
#' select_results <- select_opt(
#'   x = opt_results,
#'   smd_groups = list(
#'     c("adenoma", "controls"),
#'     c("controls", "crc_beningn"),
#'     c("crc_malignant", "controls")
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
  ## now select the desired percent_matched from the data.frame
  perc_df <- attr(x, "opt_results")

  ## select desired columns
  perc_df <- perc_df[, c("iter_ID", perc_colnames)]

  ## filter out NAs
  perc_df <- perc_df[stats::complete.cases(perc_df), ]

  ## if single column after iter_ID, then no need to compute similarity/
  ## deviance metrics - we can merge the datasets directly and proceed with
  ## selection
  if (ncol(perc_df) == 2) {
    if (colnames(perc_df)[2] != "perc_matched") {
      colnames(perc_df)[2] <- "perc_matched"
    }

    perc_df
    final_df <- merge(agg_df,
      perc_df,
      by = "iter_ID"
    )
  } else {
    ## if more than two columns calculate sum of the percentages
    # get columns after iter_ID
    metric_cols <- setdiff(names(perc_df), "iter_ID")

    # Aggregate sums by iter_ID
    perc_df$perc_matched <- rowSums(perc_df[metric_cols])

    # merge
    final_df <- merge(agg_df,
      perc_df,
      by = "iter_ID"
    )
  }

  # filter out NAs and INfs
  final_df <- final_df[apply(
    final_df, 1,
    function(row) {
      all(!is.na(row) & !is.infinite(row))
    }
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

  # define custom S3 class and attr
  class(best_rows_final) <- c("select_result", class(best_rows_final))

  # assign attribute (parameter set)
  full_df <- attr(x, "opt_results")
  remove_cols <- c( # paste0("p_", treat_names),
    "perc_matched", "smd", "smd_group"
  )

  param_df <- merge(full_df[, colnames(full_df) %nin% remove_cols],
    best_rows_final[
      ,
      colnames(best_rows_final) %nin% "perc_matched"
    ],
    by = "iter_ID"
  )

  attr(best_rows_final, "param_df") <- param_df

  # print with method
  print(best_rows_final)

  return(invisible(best_rows_final))
}

#' @export
print.select_result <- function(x, digits = 3, ...) {
  cat("\n")
  cat("Optimization Selection Summary\n")
  cat("====================\n\n")

  # Remove rows with NA in key columns
  x <- x[!is.na(x$smd_group) &
    !is.na(x$overall_stat) &
    !is.na(x$perc_matched), ]

  # Split by group
  groups <- split(x, x$smd_group)

  # Keep only groups with valid overall_stat values
  groups <- groups[sapply(groups, function(g) any(!is.na(g$overall_stat)))]

  # Sort groups by first valid overall_stat
  smd_order <- sapply(groups, function(g) min(g$overall_stat, na.rm = TRUE))
  sorted_groups <- names(sort(smd_order))

  # Build result table
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
    first_row <- rows[1, ]

    result_table <- rbind(result_table, data.frame(
      smd_group = grp,
      unique_configs = unique_configs,
      smd = round(first_row$overall_stat, digits),
      perc_matched = round(first_row$perc_matched, digits),
      stringsAsFactors = FALSE
    ))
  }

  # Table formatting
  header <- c("smd_group", "unique_configs", "smd", "perc_matched")
  col_widths <- c(12, 16, 8, 14)
  total_width <- sum(col_widths) + length(col_widths) + 1

  cat(strrep("-", total_width), "\n")
  cat("|", paste(mapply(function(name, width) {
    format(name, width = width, justify = "centre")
  }, header, col_widths), collapse = "|"), "|\n", sep = "")
  cat(strrep("-", total_width), "\n")

  for (i in seq_len(nrow(result_table))) {
    row <- result_table[i, ]
    cat("|",
      format(row$smd_group, width = col_widths[1], justify = "left"), "|",
      format(row$unique_configs, width = col_widths[2], justify = "right"), "|",
      format(row$smd, width = col_widths[3], justify = "right"), "|",
      format(row$perc_matched, width = col_widths[4], justify = "right"), "|\n",
      sep = ""
    )
  }

  cat(strrep("-", total_width), "\n\n")

  invisible(x)
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
#' - The total number of unique parameter combinations (i.e., the size of the search space).
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
    max_controls = 1) {
  ############## PARAMETER CHECK ###############################################
  ## gps_method
  ### define allowed methods
  allowed_gps_methods <- paste0("m", 1:10)
  validate_optarg(
    gps_method,
    allowed_gps_methods
  )

  ## reference
  ### define allowed treatment levels
  if (!is.null(data)) .check_df(data)

  # formula
  data_list <- .process_formula(formula, data)

  # reference processing
  treatment_var <- data_list[["treat"]]

  # args assignment to list used in calculations
  treatment_levels <- as.character(unique(treatment_var))

  # reference
  ## set default if null
  if (is.null(reference)) reference <- treatment_levels[1]

  ## process
  ref_list <- lapply(reference, function(ref) {
    .process_ref(treatment_var,
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
  ### check if numeric vector
  .chk_cond(
    !.check_vecl(caliper),
    "The argument `caliper` has to be a numeric vector."
  )

  ### check if no duplicates
  .chk_cond(
    any(duplicated(caliper)),
    "Duplicates inside the `caliper` vector are not allowed."
  )

  ### check if no negative values
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
  allowed_clusters <- 1:length(treatment_levels)

  validate_optarg(cluster,
    allowed_clusters,
    quotes = FALSE
  )

  ## replace, ties and ratio
  if ("nnm" %in% matching_method) {
    # ratio
    allowed_ratio <- 1:20
    validate_optarg(ratio,
      allowed_ratio,
      quotes = FALSE
    )

    # replace
    allowed_replace <- c(TRUE, FALSE)
    validate_optarg(
      replace,
      allowed_replace
    )

    # ties
    allowed_ties <- c(TRUE, FALSE)
    validate_optarg(
      ties,
      allowed_ties
    )
  }

  ## max and min_controls
  if ("fullopt" %in% matching_method) {
    allowed_controls <- 1:20

    validate_optarg(min_controls,
      allowed_controls,
      quotes = FALSE
    )

    validate_optarg(max_controls,
      allowed_controls,
      quotes = FALSE
    )

    ## make sure that max_controls is higher than min_controls
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
  # Collect user-specified args
  user_args <- list(
    gps_method = gps_method,
    reference = reference,
    matching_method = matching_method,
    caliper = caliper,
    order = order,
    cluster = cluster,
    replace = replace,
    ties = ties,
    ratio = ratio,
    min_controls = min_controls,
    max_controls = max_controls
  )

  # Remove anything NULL (e.g., reference not specified yet)
  user_args <- user_args[!vapply(user_args, is.null, logical(1))]

  ############# CALCULATE NUMBER OF COMBINATIONS ###############################
  ## if both methods defined, calculate separately and add
  ## otherwise, simply multiply all
  if (length(unique(user_args[["matching_method"]])) == 2) {
    total_combinations <- 0

    for (method in user_args[["matching_method"]]) {
      # Copy of the list
      method_args <- user_args

      # Exclude certain parameters based on method
      if (method == "nnm") {
        method_args[c("min_controls", "max_controls")] <- NULL
      } else if (method == "fullopt") {
        method_args[c("replace", "ties", "ratio")] <- NULL
      } else {
        warning(sprintf("Unknown matching method: '%s'. Skipping.", method))
        next
      }

      # Also fix matching_method to only this method
      method_args$matching_method <- method

      # Compute unique lengths
      unique_lengths <- vapply(
        method_args,
        function(x) length(unique(x)),
        integer(1)
      )
      total_combinations <- total_combinations + prod(unique_lengths)
    }
  } else {
    # Get lengths of unique values for each parameter
    unique_lengths <- vapply(
      user_args,
      function(x) length(unique(x)),
      integer(1)
    )

    # Calculate total combinations (cartesian product)
    total_combinations <- prod(unique_lengths)
  }

  attr(user_args, "total_combinations") <- format(total_combinations,
    scientific = FALSE
  )

  attr(user_args, "model_covs") <- colnames(data_list[["model_covs"]])

  # Set class
  class(user_args) <- "opt_args"
  return(user_args)
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
