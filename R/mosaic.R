#' Plot the distribution of categorical covariates
#'
#' @param data
#' @param y
#' @param group
#' @param facet
#' @param ncol
#' @param significance
#' @param plot.name
#' @param overwrite
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @export
mosaic <- function(data = NULL,
                   y = NULL,
                   group = NULL,
                   facet = NULL,
                   ncol = 1,
                   significance = FALSE,
                   plot.name = NULL,
                   overwrite = FALSE,
                   ...) {
  ############################ INPUT CHECKING###################################

  args_signif <- list(...)

  #--check data frame-----------------------------------------------------------
  ## convert to data.frame (to get rid of other classes passing the is.data.frame())
  class_data <- class(data)
  if ((length(class_data) > 1 && "data.frame" %in% class_data) || is.matrix(data)) {
    tryCatch(
      {
        data <- as.data.frame(data)
      },
      error = function(e) {
        chk::abort_chk("The `data` argument can not be converted to valid data.frame")
      }
    )
  }

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

  ## Check logicals
  chk::chk_all(c(overwrite, significance), chk::chk_flag)

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
    type = list("factor", "factor", "factor"),
    varname = symlist,
    MoreArgs = list(
      data = data,
      env = environment()
    )
  )

  # defining important conditions
  use_facet <- ifelse(is.null(symlist[['facet']]), FALSE, TRUE)
  use_group <- ifelse(is.null(symlist[['group']]), FALSE, TRUE)
  use_signif <- significance

  # defining levels of facet for the test
  facet_levels <- NULL
  if(use_facet) {
    facet_levels <- nunique(data[, symlist[["facet"]]])
  }

  if (facet_levels == 0 || is.null(facet_levels)) facet_levels <- 1

  # check and process the significance argument
  if (use_signif) {
    rlang::check_installed(c("rstatix", 'productplots'))

    if(!use_group) {
      chk::abort_chk('The `group` argument has to be specified if `significance = TRUE`')
    }

    # add args to the list
    args_signif[["x"]] <- data[, symlist[["group"]]]
    args_signif[["y"]] <- data[, symlist[["y"]]]

    # check groups
    if (nunique(args_signif[["x"]]) <= 1) {
      chk::abort_chk("It is impossible to compute statistical significance tests for only one group")
    }

    # perform the test for rstatix
    # matching args from ...
    args_signif <- match_add_args(
      arglist = args_signif,
      rstatix::chisq_test
    )

    # Predefine output lists
    test_results <- vector("list", facet_levels)
    res <- vector("list", facet_levels)

    # Iterate through facet levels
    for (i in seq_len(facet_levels)) {
      # Subset data for facets if there are multiple levels
      if (use_facet) {
        subset_cond <- data[, symlist[["facet"]]] == levels(data[, symlist[["facet"]]])[i]
        args_signif[c('x', 'y')] <- lapply(args_signif[c('x', 'y')],
                                           function(x) x[subset_cond])
      }

      # Run the chi-squared test and store results
      tryCatch({
        suppressWarnings({
          test_results[[i]] <- do.call(
            rstatix::chisq_test,
            args_signif
          )

          res[[i]] <- as.data.frame(rstatix::std_residuals(test_results[[i]]))
        })},
        error = function(e) {
          chk::abort_chk(sprintf(
            "There was a problem in estimating the significance levels. It is probably a bug - consider reporting to the maintainer. \n\nError message from `%s`: %s",
            "rstatix::chisq_test", conditionMessage(e)
          ), tidy = FALSE)
        }
      )

      # adding the facet var
      if (use_facet) {
        facet_name <- levels(data[, symlist[["facet"]]])[i]
        test_results[[i]] <- cbind(facet = facet_name, test_results[[i]])
        res[[i]] <- cbind(facet = facet_name, res[[i]])
      }
    }

    # process the resulting list to df
    test_results <- do.call(rbind, test_results)
    res <- do.call(rbind, res)

    # rename colnames in res
    rename_ind <- match(c('x', 'y'), colnames(res))
    colnames(res)[rename_ind] <- c(symlist[['group']],
                                   symlist[['y']])
  }

  #####################CALCULATING COORDS ######################################

  ## Define the formula (conditionally on group)
  form <- as.formula(paste("~",
                           if (use_group) paste(symlist[["group"]], '+') else NULL,
                           symlist[["y"]]))

  ## converting data to list (easier to loop on)
  data_split <- split(data,
                      if(use_facet) data[, symlist[["facet"]]] else 'no_facet')

  ## Calculate the mosaicplots coords
  #' @importFrom productplots vspine hspine
  prodcoords <- lapply(seq_along(data_split), function(i) {
    coords <- productplots::prodcalc(data_split[[i]],
                                     form,
                                     divider = productplots::mosaic(),
                                     cascade = 0,
                                     scale_max = TRUE,
                                     na.rm = FALSE
    )

    ## add facet var if necessary
    if (use_facet) {
      coords$facet <- names(data_split)[i]
    }

    coords
  })

  prodcoords <- as.data.frame(do.call(rbind, prodcoords))

  ## filter out base levels if group defined
  prodcoords <- prodcoords[complete.cases(prodcoords), ]


  ####################### PLOTTING #############################################
  ## Unique values in grouping variables (necessary to define the palette)
  pal_len <- nunique(data[, if(use_group) symlist[['group']] else symlist[['y']]])

  ## --defining the main ggplot formula------------------------------------------
  main <- ggplot2::ggplot(prodcoords)

  if(!use_group) {
    main_layers <- main +
      # defining rectangle aesthetics
      ggplot2::geom_rect(ggplot2::aes(xmin = l, xmax = r, ymin = b, ymax = t,
                                      fill = prodcoords[, symlist[["y"]]]),
                         color = 'black') +
      # customizing fill scale
      scale_fill_vecmatch(pal_len, type = "discrete") +
      # adding legend and labs
      ggplot2::guides(fill = ggplot2::guide_legend(symlist[["y"]])) +
      ggplot2::ylab(symlist[["y"]]) +
      # removing x scale
      ggplot2::theme_classic() %+replace%
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  } else {
    main <- main +
      # define ylab and xlab
      ggplot2::xlab(symlist[["y"]]) +
      # add labels on x axis
      ggplot2::ylab(symlist[["group"]]) +
      # theme
      ggplot2::theme_classic()

    if(use_signif) {
      # process test_results
      caption_text <- sprintf('\U03C7\U00B2(%i, \U1D45B = %i) = %s, \U1D45D = %s',
                              test_results$df, test_results$n,
                              format(round(test_results$statistic, digits=2), nsmall=2),
                              format(round(test_results$p, digits=3), nsmall=3))
      caption_df <- data.frame(caption = paste0('Test of independence: ',
                                               caption_text),
                               x = rep(0.6, length(caption_text)),
                               y = if(!use_facet) rep(-0.13, length(caption_text)) else rep(-0.38, length(caption_text)))

      if(use_facet) caption_df$facet <- test_results$facet

      # define merging cols
      merge_by <- if(use_facet) c('facet', symlist[['group']], symlist[['y']])
      else c(symlist[['group']], symlist[['y']])

      # join res and prodcoords to produce fill aes by standardized pearson's res
      cols_order <- colnames(prodcoords)
      prodcoords <- merge(prodcoords, res,
                          by = merge_by,
                          all.x = TRUE)

      # reorder prodcoords
      prodcoords <- prodcoords[c(cols_order, setdiff(colnames(prodcoords), cols_order))]

      # categorize standardized pearson's residuals
      sig_labels <- c('*** neg.', '** neg.', '* neg.', 'ns.',
                      '* pos.', '** pos.', '*** pos.')
      prodcoords$Freq <- as.factor(cut(prodcoords$Freq,
                             breaks = c(-Inf, -3.291, -2.576, -1.96, 1.96, 2.576, 3.291, Inf),
                             labels = sig_labels))
      fill_vals <- rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
                         "#E0F3F8", "#91BFDB", "#4575B4"))
      names(fill_vals) <- sig_labels

      # customizing the fill scale
      main_layers <- main +
        # defining rectangle aesthetics
        ggplot2::geom_rect(ggplot2::aes(xmin = l, xmax = r, ymin = b, ymax = t,
                                        fill = prodcoords[, 'Freq']),
                           color = 'black') +
        ## adding custom scale to signalize significance
        ggplot2::scale_fill_manual(name = 'Partial significance',
                                   values = fill_vals) +
        ## add caption with test results
        ggplot2::coord_cartesian(clip = 'off',
                                 ylim = c(0, 1)) +
        ggplot2::geom_text(data = caption_df,
                           ggplot2::aes(x = x, y = y, label = caption))

      if(use_facet) {
        main_layers <- main_layers +
          ggplot2::theme_classic() %+replace%
          ggplot2::theme(plot.margin = ggplot2::margin(t = 0.2, b = 1, unit = 'cm'),
                         panel.spacing = ggplot2::unit(1.5, 'cm'))
      } else {
        main_layers <- main_layers +
          ggplot2::theme_classic() %+replace%
          ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 1, 0), 'cm'))
      }

    } else {
      main_layers <- main +
        # defining rectangle aesthetics
        ggplot2::geom_rect(ggplot2::aes(xmin = l, xmax = r, ymin = b, ymax = t,
                                        fill = prodcoords[, symlist[["group"]]]),
                           color = 'black') +
        # customizing fill scale
        scale_fill_vecmatch(pal_len, type = "discrete") +
        ggplot2::guides(fill = ggplot2::guide_legend(symlist[["group"]]))
    }
  }

  # facetting if facet specified
  if (use_facet) {
    prodcoords_facet <- split(prodcoords, prodcoords$facet)

    # Create custom scales for x and y
    scales_custom_x <- if (use_group) {
      lapply(prodcoords_facet, scale_x_product)
    } else {
      vector("list", length(prodcoords_facet))  # Placeholder for x scales
    }

    scales_custom_y <- lapply(prodcoords_facet, scale_y_product)

    # Apply scale_facet to both x and y scales
    for (i in 1:facet_levels) {
      if (use_group) {
        scales_custom_x[[i]] <- scale_facet(i, scales_custom_x[[i]])
      }
      scales_custom_y[[i]] <- scale_facet(i, scales_custom_y[[i]])
    }

    # Combine scales
    scales_custom <- c(scales_custom_x, scales_custom_y)

    # facetting using custom facet_wrap
    main_layers <- main_layers +
      facet_wrap_scales(
        ~facet,
        ncol = ncol,
        scales = "free",
        scales_custom = scales_custom
      )
  } else {
    # Update main_layers based on use_group
    main_layers <- main_layers + scale_y_product(prodcoords)

    if (use_group) {
      main_layers <- main_layers + scale_x_product(prodcoords)
    }
  }
  return(main_layers)
}
