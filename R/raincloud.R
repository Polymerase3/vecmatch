#' @title Examine the imbalance of continuous covariates
#'
#' @description The `raincloud()` function allows to generate distribution plots
#'   for continuous data in an easy and uncomplicated way. The function is based
#'   on the `ggplot2` package, which must already be preinstalled Raincloud
#'   plots consist of three main elements:
#'
#' * Distribution plots, specifically  violin plots with the mean values and
#'   standard deviations of respective groups,
#' * Jittered point plots depicting the underlying distribution of the data in
#'   the rawest form,
#' * Boxplots, summarizing the most important statistics of the underlying
#'   distribution.
#'
#' @param data A non-empty `data.frame` containing at least one numeric column,
#'   as specified by the `y` argument. This argument must be provided and does
#'   not have a default value.
#' @param y A single string or unquoted symbol representing the name of a
#'   numeric column in the `data`. In the vector matching workflow, it is
#'   typically a numeric covariate that requires balancing.
#' @param group A single string or unquoted symbol representing the name of a
#'   factor or character column in `data`. In `raincloud()` plots, the groups
#'   specified by `group` argument will be distinguished by separate `fill` and
#'   `color` aesthetics. For clarity, it is recommended to plot fewer than 10
#'   groups, though there is no formal limit.
#' @param facet A single string or unquoted symbol representing the name of a
#'   variable in `data` to facet by. This argument is used in a call to
#'   [ggplot2::facet_wrap()], creating separate distribution plots for each
#'   unique group in the `facet` variable.
#' @param ncol A single integer. The value should be less than or equal to the
#'   number of unique categories in the `facet` variable. This argument is used
#'   only when `facet` is not NULL, specifying the number of columns in the
#'   [ggplot2::facet_wrap()] call. The distribution plots will be arranged into
#'   the number of columns defined by `ncol`.
#' @param significance A single string specifying the method for calculating
#'   p-values in multiple comparisons between groups defined by the `group`
#'   argument. Significant comparisons are represented by bars connecting the
#'   compared groups on the left side of the boxplots. Note that if there are
#'   many significant tests, the plot size may adjust accordingly. For available
#'   methods refer to the *Details* section. If the `significance` argument is
#'   not `NULL`, standardized mean differences (SMDs) are also calculated and
#'   displayed on the right side of the jittered point plots.
#' @param sig_label_size An integer specifying the size of the significance and
#'   SMD (standardized mean difference) labels displayed on the bars on the
#'   right side of the plot.
#' @param sig_label_color Logical flag. If `FALSE` (default), significance and
#'   SMD bars and text are displayed in the default color (black). If `TRUE`,
#'   colors are applied dynamically based on value: nonsignificant tests and SMD
#'   values below 0.10 are displayed in green, while significant tests and SMD
#'   values of 0.10 or higher are displayed in red.
#' @param smd_type A single string indicating the type of effect size to
#'   calculate and display on the left side of the jittered point plots:
#' * `mean` - Cohen's d is calculated,
#' * `median` - the Wilcoxon effect size (r) is calculated based on the Z
#'   statistic extracted from the Wilcoxon test.
#' @param limits A numeric atomic vector of length two, specifying the `y` axis
#'   limits in the distribution plots. The first element sets the minimum value,
#'   and the second sets the maximum. This vector is passed to the
#'   [ggplot2::xlim()] function to adjust the axis scale.
#' @param jitter A single numeric value between 0 and 1 that controls the amount
#'   of jitter applied to points in the [ggplot2::geom_jitter()] plots. Higher
#'   values of the `jitter` argument produce more jittered plot. It's
#'   recommended to keep this value low, as higher jitter can make the plot
#'   difficult to interpret.
#' @param alpha A single numeric value between 0 and 1 that controls the
#'   transparency of the density plots, boxplots, and jittered point plots.
#'   Lower values result in higher transparency. It is recommended to keep this
#'   value relatively high to maintain the interpretability of the plots when
#'   using the `group` argument, as excessive transparency may cause overlap
#'   between groups, making it difficult to distinguish them visually.
#' @param plot_name A string specifying a valid file name or path for the plot.
#'   If set to `NULL`, the plot is displayed to the current graphical device but
#'   not saved locally. If a valid name with `.png` or `.pdf` extension is
#'   provided, the plot is saved locally. Users can also include a subdirectory
#'   in `plot_name`. Ensure the file path follows the correct syntax for your
#'   operating system.
#' @param overwrite A logical flag (default `FALSE`) that is evaluated only if
#'   the `save.name` argument is provided. If `TRUE`, the function checks
#'   whether a plot with the same name already exists. If it does, the existing
#'   plot will be overwritten. If `FALSE` and a plot with the same name exists,
#'   an error is thrown. If no such plot exists, the plot is saved normally.
#' @param ... Additional arguments passed to the function for calculating
#'   p-values when the `significance` argument is specified. For available
#'   functions associated with different `significance` methods, please refer to
#'   the *Details* section and consult the documentation for the relevant
#'   functions in the `rstatix` package.
#'
#' @details  Available methods for the argument `significance` are:
#'* `"t_test"` - Performs a pairwise comparison using the two-sample t-test,
#'  with the default Holm adjustment for multiple comparisons. This test assumes
#'  normally distributed data and equal variances. The adjustment can be
#'  modified via the `p.adjust.method` argument. The test is implemented via
#'  [rstatix::pairwise_t_test()]
#'* `"dunn_test"` - Executes Dunn's test for pairwise comparisons following a
#'  Kruskal-Wallis test. It is a non-parametric alternative to the t-test when
#'  assumptions of normality or homogeneity of variances are violated.
#'  Implemented via [rstatix::dunn_test()].
#'* `"tukeyHSD_test"` - Uses Tukey's Honest Significant Difference (HSD) test
#'  for pairwise comparisons between group means. Suitable for comparing all
#'  pairs when the overall ANOVA is significant. The method assumes equal
#'  variance between groups and is implemented via [rstatix::tukey_hsd()].
#'* `"games_howell_test"` - A post-hoc test used after ANOVA, which does not
#'  assume equal variances or equal sample sizes. It’s particularly robust for
#'  data that violate homogeneity of variance assumptions. Implemented via
#'  [rstatix::games_howell_test()].
#'* `"wilcoxon_test"` - Performs the Wilcoxon rank-sum test (also known as the
#'  Mann-Whitney U test) for non-parametric pairwise comparisons. Useful when
#'  data are not normally distributed. Implemented via
#'  [rstatix::pairwise_wilcox_test()].
#'
#' @return A `ggplot` object representing the distribution of the `y` variable
#'   across the levels of the `group` and `facet` variables in `data`.
#'
#' @seealso [mosaic()] which summarizes the distribution of discrete data
#'
#' @examples
#' ## Example: Creating a raincloud plot for the ToothGrowth dataset.
#' ## This plot visualizes the distribution of the `len` variable by
#' ## `dose` (using different colors) and facets by `supp`. Group
#' ## differences by `dose` are calculated using a `t_test`, and standardized
#' ## mean differences (SMDs) are displayed through jittered points.
#' library(ggplot2)
#' library(ggpubr)
#'
#' p <- raincloud(ToothGrowth, len, dose, supp,
#'   significance = "t_test",
#'   jitter = 0.15, alpha = 0.4
#' )
#'
#' ## As `p` is a valid `ggplot` object, we can manipulate its
#' ## characteristics usingthe `ggplot2` or `ggpubr` packages
#' ## to create publication grade plot:
#' p <- p +
#'   theme_classic2() +
#'   theme(
#'     axis.line.y = element_blank(),
#'     axis.ticks.y = element_blank()
#'   ) +
#'   guides(fill = guide_legend("Dose [mg]")) +
#'   ylab("Length [cm]")
#'
#' p
#'
#' @export

raincloud <- function(data = NULL,
                      y = NULL,
                      group = NULL,
                      facet = NULL,
                      ncol = 1,
                      significance = NULL,
                      sig_label_size = 2L,
                      sig_label_color = FALSE,
                      smd_type = "mean",
                      limits = NULL,
                      jitter = 0.1,
                      alpha = 0.4,
                      plot_name = NULL,
                      overwrite = FALSE,
                      ...) {
  ############################ INPUT CHECKING###################################

  args_signif <- list(...)

  #--check data frame-----------------------------------------------------------
  if ("matched" %in% class(data)) {
    class(data) <- "data.frame"
  }

  ## must be an object of class data frame
  .check_df(data)

  ## with at least one numeric column
  .chk_cond(
    length(data) == 1 && !is.numeric(data[, 1]),
    "The provided data is not numeric"
  )

  #--check y, group and facet---------------------------------------------------
  ## check if the provided names are valid names + convert to characters
  symlist <- list(
    y = substitute(y),
    group = substitute(group),
    facet = substitute(facet)
  )
  symlist <- .conv_nam(symlist)

  ## check if y exists
  .chk_cond(
    is.null(symlist[[1]]),
    "The argument `y` is missing with no default!"
  )

  ## check if the names are in the provided data.frame
  nonames <- .check_name(data, symlist)
  .chk_cond(
    length(nonames) != 0,
    sprintf("The following names are not present in the
                    provided data frame: %s", word_list(add_quotes(nonames)))
  )

  ## check if limits is a numeric vector of length 2
  .chk_cond(
    !is.null(limits) && !.check_vecl(limits, leng = 2),
    "The `limits` argument should be a numeric vector of
            length 2: c(min, max)"
  )

  ## check range for jitter
  chk::chk_range(jitter, range = c(0, 1))

  ## check range for alpha
  chk::chk_range(alpha, range = c(0, 1))

  ## check logicals
  chk::chk_all(c(overwrite, sig_label_color), chk::chk_flag)

  ## check character and valid name for plot_name
  if (!is.null(plot_name)) {
    chk::chk_character(plot_name)
    .check_extension(plot_name,
      x_name = "plot_name",
      ext_vec = c(".png", ".PNG", ".pdf", ".PDF")
    )
  }

  # check if sig_label_size is integer
  if (!is.null(sig_label_size)) {
    suppressWarnings(sig_label_size <- try(as.integer(sig_label_size)))

    .chk_cond(
      is.na(sig_label_size) || inherits(sig_label_size, "try-error"),
      "`sig_label_size` must be an integer."
    )

    chk::chk_integer(sig_label_size, x_name = "sig_label_size")
  }

  # check smd_type
  chk::chk_character(smd_type)
  chk::chk_length(smd_type, length = 1)
  .chk_cond(
    smd_type %nin% c("mean", "median"),
    'The `smd_type` argument can only take one of the
            following values: "mean", "median".'
  )

  ####################### DATA PROCESSING ######################################
  # assure y is numeric and convert facet, group to factors
  mapply(.conv_data,
    type = list("numeric", "factor", "factor"),
    varname = symlist,
    MoreArgs = list(
      data = data,
      env = environment() ## allows replacing the data without reassignment
    )
  )

  ## use only complete.cases of variables in the function call
  which_use <- unlist(symlist[!vapply(symlist, is.null, logical(1L))])
  complete_sub <- stats::complete.cases(data[, colnames(data) %in% which_use])
  data <- as.data.frame(subset(data, complete_sub))

  # defining levels of facet for the test
  facet_levels <- length(unique(data[, symlist[["facet"]]]))
  facet_levels <- max(facet_levels, 1) # if facet_levels = 0

  # check and process the significance argument
  use_signif <- FALSE
  if (!is.null(significance)) {
    rlang::check_installed(c("rstatix", "ggpubr"))
    use_signif <- TRUE

    # check if more than 2 groups in the data
    .chk_cond(
      length(unique(data[, symlist[["group"]]])) <= 1,
      "It is impossible to compute statistical significance
              tests for only one group. Check your `group` argument."
    )

    # check if significance is a valid implemented method
    .chk_cond(
      !chk::vld_string(significance) || significance %nin%
        names(.sig_methods),
      sprintf(
        "The argument `significance` must be a single string
                      specifying one of the available methods: %s",
        word_list(add_quotes(names(.sig_methods)))
      )
    )

    # build formula for further calculations
    args_signif[["formula"]] <- stats::as.formula(paste0(
      symlist[["y"]],
      " ~ ",
      symlist[["group"]]
    ))

    # perform the tests for rstatix
    # define rstatix funciton used
    func_used <- switch(significance,
      "t_test" = rstatix::pairwise_t_test,
      "dunn_test" = rstatix::dunn_test,
      "tukeyHSD_test" = rstatix::tukey_hsd,
      "games_howell_test" = rstatix::games_howell_test,
      "wilcoxon_test" = rstatix::pairwise_wilcox_test,
      "sign_test" = rstatix::pairwise_sign_test
    )

    # matching args from ...
    args_signif <- match_add_args(
      arglist = args_signif,
      .sig_methods[[significance]]$args_check_fun
    )

    suppressWarnings(args_signif <- Filter(
      function(x) !all(is.na(x)),
      args_signif
    ))

    # correcting pool.sd to logical if default value for t_test
    if (!is.logical(args_signif[["pool.sd"]]) && significance == "t_test") {
      args_signif[["pool.sd"]] <- !args_signif[["paired"]]
    }

    # predefining output list
    test_results <- list()

    # fitting
    for (i in 1:facet_levels) {
      # Subsetting the data
      args_signif[["data"]] <- if (facet_levels == 1) {
        data
      } else {
        subset_cond <- data[, symlist[["facet"]]] ==
          levels(data[, symlist[["facet"]]])[i]
        data[subset_cond, ]
      }

      # modify name of data argument for tukeyHSD
      names(args_signif)[names(args_signif) == "data"] <- {
        ifelse(significance == "tukeyHSD_test", "x", "data")
      }

      ## call the rstatix func to calculate significance levels
      tryCatch(
        {
          test_results[[i]] <- do.call(
            func_used,
            args_signif
          )
        },
        error = function(e) {
          chk::abort_chk(
            strwrap(
              sprintf(
                "There was a problem in estimating the significance levels
            using %s method. It is probably a bug - consider reporting to
            the maintainer. \n
            Error message from `%s`: %s",
                significance,
                as.character(func_used),
                conditionMessage(e)
              ),
              prefix = " ", initial = ""
            ),
            tidy = FALSE
          )
        }
      )

      # calculating original add_xy_position to locate the pvalues
      test_results[[i]] <- rstatix::add_xy_position(test_results[[i]],
        fun = "max",
        stack = FALSE,
        x = symlist[["group"]],
        scales = "fixed"
      )

      # adding effsize
      names(args_signif)[names(args_signif) == "x"] <- "data"

      # dynamically choose the function based on smd_type
      effsize_function <- switch(smd_type,
        "mean" = rstatix::cohens_d,
        "median" = rstatix::wilcox_effsize
      )

      # calculate SMD using the selected function
      smd <- effsize_function(
        args_signif[["data"]],
        args_signif[["formula"]]
      )

      smd <- smd[, c("group1", "group2", "effsize")]

      colnames(smd)[colnames(smd) == "effsize"] <- "smd"

      smd$smd <- abs(round(smd$smd, 2))

      # binding the results together
      test_results[[i]] <- merge(test_results[[i]], smd,
        by = c("group1", "group2"),
        all = TRUE
      )

      # adding the facet var
      if (facet_levels > 1) {
        test_results[[i]] <- cbind(
          facet = rep(
            levels(data[, symlist[["facet"]]])[i],
            dim(test_results[[i]])[1]
          ),
          test_results[[i]]
        )
      }
    }

    # process the resulting list to df
    if (facet_levels == 1) {
      test_results <- as.data.frame(test_results)
    } else {
      test_results <- do.call(rbind, test_results)
    }

    # add color to the plot
    if (sig_label_color) {
      test_results <- cbind(test_results,
        color.pval = as.character(ifelse(test_results[, "p.adj"] < 0.05,
          "#de2d26",
          "#31a354"
        )),
        color.smd = as.character(ifelse(test_results[, "smd"] >= 0.1,
          "#de2d26",
          "#31a354"
        ))
      )
    } else {
      test_results <- cbind(test_results,
        color.pval = "#000000",
        color.smd = "#000000"
      )
    }
  }

  ####################### PLOTTING #############################################
  # define the replace function
  "%+replace%" <- ggplot2::"%+replace%"

  ## Unique values in grouping variables (necessary to define the palette)
  pal_len <- length(unique(data[, symlist[["group"]]]))
  pal_len <- max(pal_len, 1)

  ## redefining 'group' labels to include sample sizes--------------------------
  if (!is.null(symlist[["group"]])) {
    if (is.null(symlist[["facet"]])) {
      freq_table <- as.data.frame(table(data[, symlist[["group"]]]))

      for (i in seq_len(nrow(freq_table))) {
        message <- sprintf(
          "%s (n = %s)",
          freq_table[i, 1],
          freq_table[i, 2]
        )

        levels(data[, symlist[["group"]]])[
          levels(data[, symlist[["group"]]]) == freq_table[i, 1]
        ] <- message
      }
    } else {
      # define the counts
      freq_table <- as.data.frame(table(
        data[, symlist[["facet"]]],
        data[, symlist[["group"]]]
      ))

      tryCatch(
        {
          # pivoting that data.frame to wide based on sex
          freq_table <- stats::reshape(freq_table,
            timevar = "Var1",
            idvar = c("Var2"),
            direction = "wide"
          )
        },
        error = function(e) {
          chk::abort_chk(strwrap(
            "Counts could not be calculated from the provided data frame.
            Please ensure that your data frame contains observations.",
            prefix = " ", initial = ""
          ))
        }
      )


      for (i in seq_len(nrow(freq_table))) {
        message <- sprintf(
          "%s (n = %s, %s)",
          freq_table[i, 1],
          freq_table[i, 2],
          freq_table[i, 3]
        )

        levels(data[, symlist[["group"]]])[
          levels(data[, symlist[["group"]]]) == freq_table[i, 1]
        ] <- message
      }
    }
  } else {
    n_len <- length(data[, symlist[["y"]]])
  }

  #--defining height of jitter plots--------------------------------------------
  rain_height <- 0.1

  ## --defining the main ggplot formula-----------------------------------------
  colnames(data)[which(colnames(data) == symlist["facet"])] <- "facet"
  main <- ggplot2::ggplot(data, ggplot2::aes(
    x = "",
    y = data[, symlist[["y"]]]
  ))

  # defining theme for the subplots
  custom_theme <- ggplot2::theme_classic() %+replace%
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold")
    )

  ## --defining the geom_jitter
  #### CASE 1 - only y, no group, no facet======================================
  main_geom_layers <- if (pal_len == 1 || is.null(symlist[["group"]])) {
    main +
      ## --defining the geom_boxplot
      ggplot2::geom_boxplot(
        width = 0.1, alpha = alpha, show.legend = FALSE, fill = "#005b99",
        position = ggplot2::position_nudge(x = -0.22)
      ) +
      ## --defining the datapoints
      ggplot2::geom_jitter(ggplot2::aes(color = "#005b99"),
        size = 2, alpha = alpha, show.legend = TRUE,
        position = ggplot2::position_jitter(width = jitter)
      ) +
      ## halfs of the violin plots
      geom_flat_violin(
        fill = "#005b99",
        trim = FALSE, alpha = alpha,
        position = ggplot2::position_nudge(x = rain_height + 0.05)
      ) +
      ## --defining the stat_summary
      ggplot2::stat_summary(
        fun.data = mean_ci, show.legend = FALSE,
        position = ggplot2::position_nudge(x = rain_height * 3)
      ) +
      ggplot2::scale_color_identity(
        guide = "legend",
        name = symlist[["y"]],
        labels = sprintf("%s (n = %s)", symlist[["y"]], n_len)
      )
  } else {
    #### CASE 2 - y + group, no facet===========================================
    main +
      ## --defining the geom_boxplot
      ggplot2::geom_boxplot(ggplot2::aes(fill = data[, symlist[["group"]]]),
        width = 0.1, alpha = alpha, show.legend = FALSE,
        position = ggpp::position_dodgenudge(width = 0.2, x = -0.22)
      ) +
      ## --defining the datapoints
      ggplot2::geom_jitter(ggplot2::aes(color = data[, symlist[["group"]]]),
        size = 2, alpha = alpha, show.legend = FALSE,
        position = ggplot2::position_jitterdodge(
          jitter.width = jitter,
          dodge.width = 0.25
        )
      ) +
      ## halfs of the violin plots
      geom_flat_violin(ggplot2::aes(fill = data[, symlist[["group"]]]),
        trim = FALSE, alpha = alpha,
        position = ggplot2::position_nudge(x = rain_height + 0.05)
      ) +
      ## --defining the stat_summary
      ggplot2::stat_summary(ggplot2::aes(color = data[, symlist[["group"]]]),
        fun.data = mean_ci, show.legend = FALSE,
        position = ggpp::position_dodgenudge(x = rain_height * 3, width = 0.1)
      ) +
      ## define the fill lab
      ggplot2::guides(fill = ggplot2::guide_legend(symlist[["group"]])) +

      scale_color_vecmatch(n = pal_len, type = "discrete") +
      scale_fill_vecmatch(n = pal_len, type = "discrete")
  }

  #--defining the ggplot object-------------------------------------------------
  p <- main_geom_layers +
    ## defining scales
    ggplot2::scale_x_discrete(
      name = "",
      expand = c(rain_height * 3.5, 0, 0, 0.62)
    ) +
    ## flipping coordinates
    ggplot2::coord_flip() +
    ## defining theme
    ggplot2::theme_classic() %+replace%
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ## define ylabs
    ggplot2::ylab(symlist[["y"]])

  #### CASE 3 - use signif TRUE=================================================
  if (use_signif) {
    ## defining the lims and calculating the significance bars if necessary
    # test run to define the limits of the plot
    violin_test <- ggplot2::ggplot(data, ggplot2::aes(
      x = "",
      y = data[, symlist[["y"]]]
    )) +
      geom_flat_violin(
        trim = FALSE, position = ggplot2::position_nudge(x = -0.5)
      )

    # getting the xlim out of the violin plot
    x_limits <- ggplot2::ggplot_build(violin_test)
    x_limits <- as.vector(x_limits$layout$panel_scales_y[[1]]$range$range)
    range <- abs(x_limits[1] - x_limits[2])
    x_max <- max(data[, symlist[["y"]]]) + range * 0.03
    x_limits[2] <- x_limits[2] + range * 0.06

    # overwrite x_limits if limits defined
    x_limits <- limits %||% x_limits

    # getting the nudged xlim
    p_build <- ggplot2::ggplot_build(p)
    x_boxplot <- p_build$data[[1]][, c("group", "x")]
    x_jitter <- p_build$data[[2]][, c("group", "x")]

    # recalculating the y.positions
    # calculating approximate range, maximum y and y.positions for stat_pvalue
    y_pos_seq <- list()

    # looping along the levels of fcaet to calulacte yposition
    for (i in 1:facet_levels) {
      #### CASE 3.1 - use facet=================================================
      if (facet_levels > 1) {
        test_results_sub <- test_results[
          test_results$facet == unique(levels(data[, "facet"]))[i],
        ]
      } else {
        #### CASE 3.1 - no facet================================================
        test_results_sub <- test_results
      }

      number_comp <- dim(test_results_sub)[1]
      y_pos_seq[[i]] <- seq(
        from = x_max,
        to = x_limits[2],
        length.out = number_comp
      )
    }

    # overwriting y.position
    test_results$y.position <- unlist(y_pos_seq)

    # overwriting xmin and xmax
    test_results$xmin_box <- x_boxplot$x[match(
      test_results$xmin,
      x_boxplot$group
    )]
    test_results$xmax_box <- x_boxplot$x[match(
      test_results$xmax,
      x_boxplot$group
    )]
    test_results$xmin_jit <- x_jitter$x[match(
      test_results$xmin,
      x_jitter$group
    )]
    test_results$xmax_jit <- x_jitter$x[match(
      test_results$xmax,
      x_jitter$group
    )]

    ## adding stat_pvalue
    p <- p +
      ## adding custom pvalues
      ggpubr::stat_pvalue_manual(test_results,
        y.position = "y.position",
        xmin = "xmin_box",
        xmax = "xmax_box",
        label = "p.adj.signif",
        coord.flip = TRUE,
        tip.length = 0.01,
        size = sig_label_size,
        color = rep(test_results[, "color.pval"], each = 3)
      ) +
      ## adding smds
      ggpubr::stat_pvalue_manual(test_results,
        y.position = "y.position",
        xmin = "xmin_jit",
        xmax = "xmax_jit",
        label = "smd",
        coord.flip = TRUE,
        tip.length = 0.01,
        size = sig_label_size,
        color = rep(test_results[, "color.smd"], each = 3)
      )
  }

  ## add theme
  p <- p +
    custom_theme

  # define limits
  if (!is.null(limits)) {
    p <- p +
      ggplot2::ylim(limits)
  }

  #### CASE 4 - add facet=======================================================
  if (!is.null(symlist[["facet"]])) {
    ## custom labeller function to include the number of obs--------------------
    freq_table <- as.data.frame(table(data[, "facet"]))

    labeller_vec <- c()
    for (i in seq_len(nrow(freq_table))) {
      message <- sprintf(
        "%s (n = %s)",
        freq_table[i, 1],
        freq_table[i, 2]
      )

      ## defining named variable
      labeller_vec[i] <- message
      names(labeller_vec)[i] <- as.character(freq_table[i, 1])
    }

    p <- p +
      ggplot2::facet_wrap(. ~ facet,
        ncol = ncol,
        labeller = ggplot2::as_labeller(labeller_vec)
      )
  }

  ## Saving the plot
  if (!is.null(plot_name)) {
    fexist <- file.exists(plot_name)

    if (overwrite || (!fexist && !overwrite)) {
      suppressMessages(ggplot2::ggsave(plot_name,
        plot = p, dpi = 300, create.dir = TRUE
      ))
    } else if (fexist && !overwrite) {
      chk::abort_chk(strwrap(
        "The file specified in the plot_name argument already exists.
        Set overwrite = TRUE to replace the existing file.",
        prefix = " ", initial = ""
      ))
    }
  }

  ## Returning a ggplot object
  return(p)
}
