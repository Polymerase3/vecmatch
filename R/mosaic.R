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
#'
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
  if((length(class_data) > 1 && 'data.frame' %in% class_data) || is.matrix(data)) {
    tryCatch({
      data <- as.data.frame(data)},
      error = function(e) {
        chk::abort_chk('The `data` argument can not be converted to valid data.frame')
      })
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

  # defining levels of facet for the test
  facet_levels <- length(unique(data[, symlist[['facet']]]))
  if(facet_levels == 0 || is.null(facet_levels)) facet_levels <- 1

  # check and process the significance argument
  if (significance) {
    rlang::check_installed(c("rstatix", "ggpubr"))

    # check groups
    if (length(unique(data[, symlist[["group"]]])) <= 1) {
      chk::abort_chk("It is impossible to compute statistical significance tests for only one group")
    }

    # add args
    args_signif[["x"]] <- data[, symlist[["group"]]]
    args_signif[['y']] <- data[, symlist[['y']]]

    # perform the test for rstatix
      # matching args from ...
      args_signif <- match_add_args(
        arglist = args_signif,
        rstatix::chisq_test
      )

      #predefining output list
      test_results <- list()
      res <- list()

      #fitting
      for (i in 1:facet_levels) {

        # subsetting the data
        if(facet_levels > 1) {
          subset_cond <- data[, symlist[['facet']]] == levels(data[, symlist[['facet']]])[i]
          args_signif[["y"]] <- args_signif$y[subset_cond]
          args_signif[['x']] <- args_signif$x[subset_cond]
        }

        ## call the rstatix func
        tryCatch(
          {
            test_results[[i]] <- do.call(
              rstatix::chisq_test,
              args_signif
            )

            res[[i]] <- rstatix::std_residuals(test_results[[i]])

          },
          error = function(e) {
            chk::abort_chk(sprintf(
              'There was a problem in estimating the significance levels using %s method. It is probably a bug - consider reporting to the maintainer. \n
            Error message from `%s`: %s',
              significance, 'rstatix::chisq_test', conditionMessage(e)
            ), tidy = FALSE)
          }
        )

        # adding the facet var
        if(facet_levels > 1) {
          test_results[[i]] <- cbind(facet = levels(data[, symlist[['facet']]])[i],
                                     test_results[[i]])
        }
      }

      # process the resulting list to df
      if(facet_levels == 1) {
        test_results <- as.data.frame(test_results)
        res <- as.data.frame(res)
      } else {
        test_results <- do.call(rbind, test_results)
        res <- lapply(res, as.data.frame)
        names(res) <- levels(data[, symlist[['facet']]])
        ##
      }
  }

  ####################### PLOTTING #############################################
  # define the replace function
  "%+replace%" <- ggplot2::"%+replace%"

  ## Unique values in grouping variables (necessary to define the palette)
  pal_len <- length(unique(data[, symlist[["group"]]]))
  if (is.null(symlist[["group"]])) {
    pal_len <- length(unique(data[, symlist[["y"]]]))
  }


  ## Define the formula
  if(!is.null(symlist[['group']])) {
    form <- as.formula(paste('~ ', symlist[['group']],  ' + ', symlist[['y']]))
  } else {
    form <- as.formula(paste('~ ', symlist[['y']]))
  }

  ## converting data to list (easier to loop on)
  if(facet_levels == 1) {
    data_split <- list(no_facet = data)
  } else {
    data_split <- split(data, data[symlist[['facet']]])
  }

  ## Define vars for iteration
  prodcoords <- list()
  i = 1

  ## Calculate the mosaicplots coords
  #' @importFrom productplots vspine hspine
  for (levels in data_split) {
    prodcoords[[i]] <- productplots::prodcalc(levels,
                                              form,
                                              divider = productplots::mosaic(),
                                              cascade = 0,
                                              scale_max = TRUE,
                                              na.rm = FALSE)
    ## add facet var if necessary
    if(facet_levels > 1) {
      prodcoords[[i]]$facet <- names(data_split)[i]
      i = i + 1
    }
  }

  ## convert output
  if(facet_levels == 1) {
    prodcoords <- as.data.frame(prodcoords)
  } else {
    prodcoords <- as.data.frame(do.call(rbind, prodcoords))
  }

  ## filter out base levels if group defined
  if(!is.null(symlist[['group']])) prodcoords <- prodcoords[prodcoords$level != 1, ]

  ## --defining the main ggplot formula------------------------------------------
  main <- ggplot2::ggplot(prodcoords) +
    # defining deafult rectangle aesthetics
    ggplot2::aes(xmin = l, xmax = r, ymin = b, ymax = t)

  ## expand to 'group' argument
  if(!is.null(symlist[['group']])) {
    main <- main +
      # adding fill aes
      ggplot2::aes(fill = prodcoords[, symlist[['group']]]) +
      # customizing fill scale
      scale_fill_vecmatch(pal_len, type = 'discrete') +
      # # add labels on x axis
      scale_x_product(prodcoords) +
      # ## define ylab and xlab and fill legend
      ggplot2::xlab(symlist[['y']]) +
      ggplot2::ylab(symlist[['group']]) +
      ggplot2::guides(fill = ggplot2::guide_legend(symlist[['group']])) +
      ggplot2::theme_classic()
  } else {
    main <- main +
      # adding fill aes
      ggplot2::aes(fill = prodcoords[, symlist[['y']]]) +
      # customizing fill scale
      scale_fill_vecmatch(pal_len, type = 'discrete') +
      # removing x scale
      ggplot2::theme_classic() %+replace%
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.line.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      # adding legend and labs
      ggplot2::guides(fill = ggplot2::guide_legend(symlist[['y']])) +
      ggplot2::ylab(symlist[['y']])
  }

  ## expand to 'facet' argument

  main <- main +
    scale_y_product(prodcoords) +
    ggplot2::geom_rect()

  if(facet_levels > 1) {
    main <- main +
      ggplot2::facet_wrap(~ facet, ncol = 1)
  }
  print(main)
}

scale_x_product <- function(coords) {
  # split data based on x axis var (always first in the productplot output)
  coords <- split(coords, coords[, 1])

  # predefining vars
  label_position <- list()
  i = 1

  # calculating label positions and defining names
  for(subs in coords) {
    coords_sub <- as.data.frame(subset(subs, subs $b == 0))
    label_position$pos[i] <- (coords_sub$l + coords_sub$r)/2
    label_position$name[i] <- names(coords)[i]
    i = i + 1
  }

  # adding scale to plot
  ggplot2::scale_x_continuous(breaks = label_position$pos,
                              labels = label_position$name)
}

scale_y_product <- function(coords) {
  # subset data (only lefts)
  coords_sub <- as.data.frame(subset(coords, coords$l == 0))

  if(!is.null(coords_sub$facet)) {
    coords_sub <- subset(coords_sub, coords_sub$facet == unique(coords_sub$facet)[1])
    coords_sub <- coords_sub[, -ncol(coords_sub)]
  }

  print(coords_sub)
  # predefining vars
  label_position <- list()
  single_aes <- unique(coords_sub$level) == 1

  # calculating the positions and labels
  for(i in 1:nrow(coords_sub)) {
    sub_data <- coords_sub[i, ]
    label_position$pos[i] <- (sub_data$b + sub_data$t) / 2

    if(single_aes) {
      label_position$name[i] <- as.character(sub_data[, 1])
    } else {
      label_position$name[i] <- as.character(sub_data[, ncol(sub_data)])
    }
  }

  # adding scale to plot
  ggplot2::scale_y_continuous(breaks = label_position$pos,
                              labels = label_position$name)
}
