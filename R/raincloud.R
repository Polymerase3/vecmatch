#' @title Examine the imbalance of continuous covariates
#'
#' @description The `raincloud()` function allows to generate distribution plots
#'   for continuous data in an easy and uncomplicated way.
#'
#'   Raincloud plots consist of three main elements: Distribution plots For
#'   example density plots, histograms or violin plots wit the mean values of
#'   respective groups Jittered point plots Depicting the underlying
#'   distribution of the data in the rawest form Boxplots The function is based
#'   on the `ggplot2` package, which must already be preinstalled
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
#'   methods refer to the Details section.
#' @param smd A logical flag, defaulting to `FALSE`. Specifies whether to
#'   display the pairwise standardized mean differences (SMD) between the groups
#'   defined by the `group` argument. SMD values are shown on the left side of
#'   the boxplots, similarly to the significance levels from the `significance`
#'   argument. Note that only one of these two arguments, `smd` or
#'   `significance`, can be set to `TRUE` in a single function call.
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
#' @param save.name A string specifying a valid file name or path for the plot.
#'   If set to `NULL`, the plot is displayed to the current graphical device but
#'   not saved locally. If a valid name with `.png` or `.pdf` extension is
#'   provided, the plot is saved locally. Users can also include a subdirectory
#'   in `save.name`, but the directory must be created manually. Ensure the file
#'   path follows the correct syntax for your operating system.
#' @param overwrite A logical flag (default `FALSE`) that is evaluated only if
#'   the `save.name` argument is provided. If `TRUE`, the function checks
#'   whether a plot with the same name already exists. If it does, the existing
#'   plot will be overwritten. If `FALSE` and a plot with the same name exists,
#'   an error is thrown. If no such plot exists, the plot is saved normally.
#'
#' @details  Available methods for the argument `significance` are:
#'  * `"t_test"` - Performs a pairwise comparison using the two-sample t-test, with the default Holm adjustment for multiple comparisons. This test assumes normally distributed data and equal variances. The adjustment can be modified via the `p.adjust.method` argument. The test is implemented via [rstatix::pairwise_t_test()]
#'  * `"dunn_test"` - Executes Dunn's test for pairwise comparisons following a Kruskal-Wallis test. It is a non-parametric alternative to the t-test when assumptions of normality are violated. Implemented via [rstatix::dunn_test()].
#'  * `"tukeyHSD_test"` - Uses Tukey's Honest Significant Difference (HSD) test for pairwise comparisons between group means. Suitable for comparing all pairs when the overall ANOVA is significant. The method assumes equal variance between groups and is implemented via [rstatix::tukey_hsd()].
#'  * `"games_howell_test"` - A post-hoc test used after ANOVA, which does not assume equal variances or equal sample sizes. It’s particularly robust for data that violate homogeneity of variance assumptions. Implemented via [rstatix::games_howell_test()].
#'  * `"wilcoxon_test"` - Performs the Wilcoxon rank-sum test (also known as the Mann-Whitney U test) for non-parametric pairwise comparisons. Useful when data are not normally distributed. Implemented via [rstatix::pairwise_wilcox_test()].
#'  * `"sign_test"` - Conducts the sign test for matched pairs, a non-parametric alternative to the paired t-test when the data don’t meet parametric assumptions. Implemented via [rstatix::pairwise_sign_test()].
#'  * `"scheffe_test"` - This test performs Scheffé’s method, a post-hoc analysis used after ANOVA. Scheffé’s method is conservative and useful for making complex comparisons. Implemented via [multcomp::glht()].
#'  * `"LSD_test"` - The Least Significant Difference (LSD) test compares group means with no adjustment for multiple comparisons, making it more powerful but more prone to Type I errors. Implemented via [multcomp::glht()].
#'  * `"sidak_test"` - Adjusts p-values for multiple comparisons using the Šidák correction, which is more powerful than the Bonferroni correction but still conservative. Implemented via [multcomp::glht()].
#'  * `"hochberg_test"` - Uses Hochberg's method for controlling the family-wise error rate in multiple comparisons. It is a step-up procedure and more powerful than Bonferroni. Implemented via [multcomp::glht()].
#'
#' @returns
#' @export
#'
#' @examples
raincloud <- function(data = NULL,
                      y = NULL,
                      group = NULL,
                      facet = NULL,
                      ncol = 1,
                      significance = NULL, ## not functional
                      smd = FALSE, # actually wywalić + merge z significance
                      limits = NULL, # not functional
                      jitter = 0.1,
                      alpha = 0.4,
                      plot.name = NULL,
                      overwrite = FALSE,
                      ...) { ## not functional
  ############################ INPUT CHECKING###################################

  args_signif <- list(...)

  #--check data frame-----------------------------------------------------------
  ## Must be an object of class data frame with at least one numeric column
  .check_df(data)

  if (length(data) == 1 && !is.numeric(data[, 1])) {
    chk::abort_chk("The provided data is not numeric")
  }

  #--check y, group and facet---------------------------------------------------

  ## Check if the provided names are valid names + convert to
  symlist <- list(
    y = substitute(y),
    group = substitute(group),
    facet = substitute(facet)
  )
  symlist <- .conv_nam(symlist)

  ## Check if y exists
  if (is.null(symlist[[1]])) chk::abort_chk("The argument `y` is missing
                                           with no default!")

  ## Check if there are in the dataframe
  nonames <- .check_name(data, symlist)

  if (length(nonames) != 0) {
    chk::abort_chk(paste0(
      "The following colnames are not in the
                          provided data frame: ",
      paste(nonames, collapse = ", ")
    ))
  }

  ## Check if limits is a numeric vector of length 2
  if (!is.null(limits)) {
    if (!.check_vecl(limits, leng = 2)) {
      chk::abort_chk("The `limits` argument should be a numeric vector of length
                     2: c(`min`, `max`)")
    }
  }

  ## Check range for jitter
  if (!is.null(jitter)) chk::chk_range(jitter, range = c(0, 1))

  ## Check range for alpha
  chk::chk_range(alpha, range = c(0, 1))

  ## Check logicals
  chk::chk_all(c(smd, overwrite), chk::chk_flag)

  ## Check character and valid name for plot.name
  if (!is.null(plot.name)) {
    chk::chk_character(plot.name)
    .check_extension(plot.name,
      x_name = "plot.name",
      ext_vec = c(".png", ".PNG", ".pdf", ".PDF")
    )
  }

  ####################### DATA PROCESSING ######################################
  # assure y is numeric and convert facet, group to factors
  mapply(.conv_data,
    type = list("numeric", "factor", "factor"),
    varname = symlist,
    MoreArgs = list(
      data = data,
      env = environment()
    )
  )

  # check and process the significance argument
  use.signif <- FALSE
  if (!is.null(significance)) {
    rlang::check_installed(c("rstatix", "multcomp", "ggpubr"))
    use.signif <- TRUE

    # check groups
    if (length(unique(data[, symlist[["group"]]])) <= 1) {
      chk::abort_chk("It is impossible to compute statistical significance tests for only one group")
    }

    # defining possible methods
    methods <- list(
      "t_test" = list(
        method_name = "t_test",
        package_used = "rstatix",
        args_check_fun = list(
          rstatix::pairwise_t_test
        )
      ),
      "dunn_test" = list(
        method_name = "dunn_test",
        package_used = "rstatix",
        args_check_fun = list(
          rstatix::dunn_test
        )
      ),
      "tukeyHSD_test" = list(
        method_name = "tukeyHSD",
        package_used = "rstatix",
        args_check_fun = list(
          utils::getS3method("tukey_hsd", "data.frame", envir = asNamespace("rstatix"))
        )
      ),
      "games_howell_test" = list(
        method_name = "games_howell_test",
        package_used = "rstatix",
        args_check_fun = list(
          rstatix::games_howell_test
        )
      ),
      "wilcoxon_test" = list(
        method_name = "wilcoxon_test",
        package_used = "rstatix",
        args_check_fun = list(
          rstatix::pairwise_wilcox_test
        )
      ),
      "sign_test" = list(
        method_name = "sign_test",
        package_used = "rstatix",
        args_check_fun = list(
          rstatix::pairwise_sign_test,
          rstatix::sign_test
        )
      ),
      "scheffe_test" = list(
        method_name = "scheffe_test",
        package_used = "multcomp",
        args_check_fun = list(
          multcomp::glht
        )
      ),
      "LSD_test" = list(
        method_name = "LSD_test",
        package = "multcomp",
        args_check_fun = list(
          multcomp::glht
        )
      ),
      "sidak_test" = list(
        method_name = "sidak_test",
        package = "multcomp",
        args_check_fun = list(
          multcomp::glht
        )
      ),
      "hochberg_test" = list(
        method_name = "hochberg_test",
        package = "multcomp",
        args_check_fun = list(
          multcomp::glht
        )
      )
    )

    # check significance
    if (!chk::vld_string(significance) || significance %nin% names(methods)) {
      chk::abort_chk(sprintf(
        "The argument `significance` has to be a single R string, describing one of the available methods: %s",
        word_list(add_quotes(methods))
      ))
    }

    # build formula
    args_signif[["formula"]] <- as.formula(paste0(symlist[["y"]], " ~ ", symlist[["group"]]))

    # args data


    # perform the tests
    if (methods[[significance]]$package_used == "rstatix") {
      # modify name of data argument for tukeyHSD
      if (significance == "tukeyHSD_test") {
        names(args_signif)[which(names(args_signif) == "data")] <- "x"
      }

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
        methods[[significance]]$args_check_fun
      )

      # correcting pool.sd to logical if default value for t_test
      if (!is.logical(args_signif[["pool.sd"]]) && significance == "t_test") {
        args_signif[["pool.sd"]] <- !args_signif[["paired"]]
      }

      # defining levels of facet for the test
      facet_levels <- length(unique(data[, symlist[['facet']]]))
      if(facet_levels == 0 || is.null(facet_levels)) facet_levels <- 1

      #predefining output list
      test_results <- list()

      #fitting
      for (i in 1:facet_levels) {

        # subsetting the data
        if(facet_levels == 1) {
          args_signif[["data"]] <- data
        } else {
          subset_cond <- data[, symlist[['facet']]] == levels(data[, symlist[['facet']]])[i]
          args_signif[["data"]] <- data[subset_cond, ]
        }

        ## call the rstatix func
        tryCatch(
          {
            test_results[[i]] <- do.call(
              func_used,
              args_signif
            )
          },
          error = function(e) {
            chk::abort_chk(sprintf(
              'There was a problem in estimating the significance levels using %s method. It is probably a bug - consider reporting to the maintainer. \n
            Error message from `%s`: %s',
              significance, as.character(func_used), conditionMessage(e)
            ), tidy = FALSE)
          }
        )

        # calculating original add_xy_position
        test_results[[i]] <- rstatix::add_xy_position(test_results[[i]],
                                                 fun = "max",
                                                 stack = FALSE, x = symlist[["group"]],
                                                 scales = "fixed")

        # adding effsize
        if (significance == "tukeyHSD_test") {
          smd <- rstatix::cohens_d(
            args_signif[["x"]],
            args_signif[["formula"]]
          )
        } else {
          smd <- rstatix::cohens_d(
            args_signif[["data"]],
            args_signif[["formula"]]
          )
        }

        # binding the results
        test_results[[i]] <- cbind(test_results[[i]], smd = abs(round(smd$effsize, 2)))

        # adding the facet var
        if(facet_levels > 1) {
          test_results[[i]] <- cbind(facet = rep(levels(data[, symlist[['facet']]])[i],
                                                 dim(test_results[[i]])[1]),
                                     test_results[[i]])
        }
      }

      # process the resulting list to df
      if(facet_levels == 1) {
        test_results <- as.data.frame(test_results)
      } else {
        test_results <- do.call(rbind, test_results)
      }

    # deal with the multcomp package and method defined in it
    } else if (attr(significance, "package") == "multcomp") {

    }
    ## calculate table
  }

  ####################### PLOTTING #############################################
  # define the replace function
  "%+replace%" <- ggplot2::"%+replace%"

  ## Unique values in grouping variables (necessary to define the palette)
  pal_len <- length(unique(data[, symlist[["group"]]]))

  if (use.signif) {
    # define the main
    main <- ggplot2::ggplot(data, ggplot2::aes(
      x = data[, symlist[["group"]]],
      y = data[, symlist[["y"]]]
    ))

    # defining theme for the subplots
    custom_theme <- ggplot2::theme_classic() %+replace%
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(face = "bold")
      )

    # test run to define the limits of the plot
    violin_test <- ggplot2::ggplot(data, ggplot2::aes(
      x = "",
      y = data[, symlist[["y"]]],
      fill = data[, symlist[["group"]]]
    )) +
      geom_flat_violin(
        trim = FALSE, alpha = alpha, show.legend = TRUE,
        position = ggplot2::position_nudge(x = -0.5)
      )

    # getting the xlim
    x_limits <- ggplot2::ggplot_build(violin_test)
    x_limits <- as.vector(x_limits$layout$panel_scales_y[[1]]$range$range)
    range <- abs(x_limits[1] - x_limits[2])
    x_max <- x_limits[2] + range * 0.03
    x_limits[2] <- x_limits[2] + range * 0.1

    # recalculating the y.positions
    # calculating approximate range, maximum y and y.positions for stat_pvalue
    if(facet_levels == 1) {
      number_comp <- dim(test_results)[1]
      print(number_comp)
      y.pos.seq <- list(seq(
        from = x_max,
        to = x_limits[2],
        length.out = number_comp
      ))
    } else {
      y.pos.seq <- list()
      # looping along the levels of fcaet to calulacte yposition
      for(i in 1:facet_levels) {
        test_results_sub <- test_results[test_results$facet == unique(levels(data[, symlist[['facet']]]))[i], ]
        number_comp <- dim(test_results_sub)[1]
        y.pos.seq[[i]] <- seq(
          from = x_max,
          to = x_limits[2],
          length.out = number_comp
        )
      }
    }

    # overwriting y.position
    test_results$y.position <- unlist(y.pos.seq)

    # generating violins
    violin_plot <- ggplot2::ggplot(data, ggplot2::aes(
      x = "",
      y = data[, symlist[["y"]]],
      fill = data[, symlist[["group"]]]
    )) +
      geom_flat_violin(
        trim = FALSE, alpha = alpha, show.legend = FALSE,
        position = ggplot2::position_nudge(x = -0.5)
      ) +
      scale_fill_vecmatch(n = pal_len, type = "discrete") +
      scale_color_vecmatch(n = pal_len, type = "discrete") +
      ## defining scales

      ggplot2::scale_x_discrete(
        name = "",
        expand = c(0, 0)
      ) +
      ## flipping coordinates
      ggplot2::coord_flip() +
      # defining theme
      custom_theme %+replace%
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      ) +
      ## --defining the stat_summary
      ggplot2::stat_summary(ggplot2::aes(color = data[, symlist[["group"]]]),
        fun.data = ggplot2::mean_cl_normal, show.legend = FALSE,
        position = ggpp::position_dodgenudge(x = -0.3, width = 0.1)
      ) +
      ggplot2::ylim(x_limits)

    # generating boxplot
    box_plot <- main +
      ggplot2::geom_boxplot(ggplot2::aes(fill = data[, symlist[["group"]]]),
        width = 0.5, alpha = alpha, show.legend = FALSE,
        position = ggplot2::position_dodge()
      ) +
      scale_fill_vecmatch(n = pal_len, type = "discrete") +
      ## flipping coordinates
      ggplot2::coord_flip() +
      custom_theme +
      ggplot2::ylim(x_limits) +
      ggplot2::ylab(symlist[["y"]])

    ## generating jitters
    jitter_plot <- main +
      ggplot2::geom_jitter(ggplot2::aes(color = data[, symlist[["group"]]]),
        size = 2, alpha = alpha, show.legend = TRUE,
        position = ggplot2::position_jitter(width = jitter)
      ) +
      scale_color_vecmatch(n = pal_len, type = "discrete") +
      ## flipping coordinates
      ggplot2::coord_flip() +
      custom_theme %+replace%
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::ylim(x_limits) +
      ggplot2::labs(color = symlist[["group"]])

    if(facet_levels == 1) {
      box_plot <- box_plot +
        ## adding custom pvalues
        ggpubr::stat_pvalue_manual(test_results,
                                   label = "p.adj.signif", # Use p.signif to display p-values as asterisks (optional)
                                   coord.flip = TRUE,
                                   tip.length = 0.01
        )

      jitter_plot <- jitter_plot +
        ggpubr::stat_pvalue_manual(test_results,
                                   label = "smd", # Use p.signif to display p-values as asterisks (optional)
                                   coord.flip = TRUE,
                                   tip.length = 0.01
        )
    }

    p <- ggpubr::ggarrange(violin_plot, jitter_plot, box_plot,
      ncol = 1,
      heights = c(2, 1, 1),
      common.legend = TRUE,
      legend = "right"
    )
  } else {
    #--defining necessary variables-----------------------------------------------
    rain_height <- 0.1

    ## --defining the main ggplot formula------------------------------------------
    main <- ggplot2::ggplot(data, ggplot2::aes(x = "", y = data[, symlist[["y"]]]))

    if (!is.null(symlist[["group"]])) {
      main <- main +
        ggplot2::aes(
          fill = data[, symlist[["group"]]],
          color = data[, symlist[["group"]]]
        )
    }

    if (is.null(symlist[["group"]]) || pal_len == 0) pal_len <- 1

    ## --defining the geom_jitter
    main_geom_layers <- if (pal_len == 1 || is.null(symlist[["group"]])) {
      main +
        ## --defining the datapoints
        ggplot2::geom_jitter(
          size = 2, alpha = alpha, show.legend = FALSE,
          position = ggplot2::position_jitter(width = jitter)
        ) +

        ## --defining the geom_boxplot
        ggplot2::geom_boxplot(
          width = 0.1, alpha = alpha, show.legend = FALSE,
          position = ggplot2::position_nudge(x = -0.22)
        ) +

        ## --defining the stat_summary
        ggplot2::stat_summary(
          fun.data = ggplot2::mean_cl_normal, show.legend = FALSE,
          position = ggplot2::position_nudge(x = rain_height * 3)
        ) +

        ## halfs of the violin plots
        geom_flat_violin(
          trim = FALSE, alpha = alpha,
          position = ggplot2::position_nudge(x = rain_height + 0.05)
        )
    } else {
      main +
        ## --defining the geom_boxplot
        ggplot2::geom_boxplot(
          width = 0.1, alpha = alpha, show.legend = FALSE,
          position = ggpp::position_dodgenudge(width = 0.2, x = -0.22)
        ) +
        ## --defining the datapoints
        ggplot2::geom_jitter(
          size = 2, alpha = alpha, show.legend = FALSE,
          position = ggplot2::position_jitterdodge(
            jitter.width = jitter,
            dodge.width = 0.25
          )
        ) +

        ## --defining the stat_summary
        ggplot2::stat_summary(
          fun.data = ggplot2::mean_cl_normal, show.legend = FALSE,
          position = ggpp::position_dodgenudge(x = rain_height * 3, width = 0.1)
        ) +

        ## halfs of the violin plots
        geom_flat_violin(
          trim = FALSE, alpha = alpha,
          position = ggplot2::position_nudge(x = rain_height + 0.05)
        )
    }

    #--defining the ggplot object-------------------------------------------------
    p <- main_geom_layers +
      ## defining scales
      ggplot2::scale_x_discrete(
        name = "",
        expand = c(rain_height * 3.5, 0, 0, 0.62)
      ) +
      scale_color_vecmatch(n = pal_len, type = "discrete") +
      scale_fill_vecmatch(n = pal_len, type = "discrete") +
      ## flipping coordinates
      ggplot2::coord_flip() +
      ## defining theme
      ggplot2::theme_classic() %+replace%
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(face = "bold")
      )
  }

  #--add facet if not NULL------------------------------------------------------
  if (!is.null(symlist[["facet"]])) {
    print(data[, symlist[["facet"]]])
    p <- p +
      ggplot2::facet_wrap(~ data[, symlist[["facet"]]], ncol = ncol)
  }

  #--save if specified
  ## Saving the plot
  # if (save == TRUE) {
  #   suppressMessages(ggplot2::ggsave(plot.name,
  #     plot = p, dpi = 300, create.dir = TRUE
  #   ))
  # }

  ## Returning a ggplot object
  return(p)
}
