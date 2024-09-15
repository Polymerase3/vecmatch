#' Title
#'
#' @return
#' @export
#'
#' @examples
raincloud <- function(data = NULL,
                      y = NULL,
                      group = NULL,
                      facet = NULL,
                      ncol = 1,
                      significance = FALSE,
                      limits = NULL,
                      jitter = 0.1,
                      alpha = 0.4,
                      save = NULL) {
  ############################ INPUT CHECKING###################################
  #--check data frame-----------------------------------------------------------
  ## Must be an object of class data frame with at least one numeric column
  if (is.null(data) || !inherits(data, "data.frame")) {
    chk::abort_chk("Argument `data` must be an object of class `data.frame`")
  }

  if (length(data) == 0) {
    chk::abort_chk("The provided data frame is empty")
  }

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

  ## Check if significance is logical
  chk::chk_logical(significance)

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

  ## Check logical for save
  chk::chk_logical(save)

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

  ####################### PLOTTING #############################################
  ##--defining the main ggplot formula------------------------------------------
  main <- paste0(
    "ggplot2::ggplot(data, ggplot2::aes(x = '', y = y",
    ifelse(is.null(symlist[["group"]]), "))", paste0(
      ", fill = ",
      symlist[["group"]], ", color = ", symlist[["group"]],
      "))"
    ))
  )

  #--defining necessary variables-----------------------------------------------
  "%+replace%" <- ggplot2::"%+replace%"
  rain_height <- 0.1

  ## Unique values in grouping variables (necessary to define the palette)
  pal_len <- length(unique(data[, symlist[["group"]]]))
  if (is.null(symlist[["group"]]) || pal_len == 0) pal_len <- 1

  #--defining the ggplot object-------------------------------------------------
  p <- eval(parse(text = main)) +
    ## halfs of the violin plots
    geom_flat_violin(
      trim = FALSE, alpha = alpha,
      position = ggplot2::position_nudge(x = rain_height + 0.05)
    ) +
    ## datapoints
    ggplot2::geom_jitter(
      size = 2, alpha = alpha, show.legend = FALSE,
      position = ggplot2::position_jitterdodge(
        jitter.width = jitter,
        dodge.width = 0.25
      )
    ) +
    ## boxplot
    ggplot2::geom_boxplot(
      width = 0.1, alpha = alpha, show.legend = FALSE,
      position = ggpp::position_dodgenudge(width = 0.2, x = -0.22)
    ) +
    ## stat_summary
    ggplot2::stat_summary(
      fun.data = ggplot2::mean_cl_normal, show.legend = FALSE,
      position = ggpp::position_dodgenudge(
        width = 0.1,
        x = rain_height * 3
      )
    ) +
    ## defining scales
    ggplot2::scale_x_discrete(name = "",
                              expand = c(rain_height * 3.5, 0, 0, 0.62)) +
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

  #--add facet if not NULL------------------------------------------------------
  if(!is.null(symlist[['facet']])) {
    p <- p +
      ggplot2::facet_wrap(~data[, symlist[['facet']]], ncol = ncol)
  }

  #--save if specified
  ## Saving the plot
  if (save == TRUE) {
    ggplot2::ggsave('plots/raincloud.png',
           plot = p, dpi = 300, create.dir = TRUE
    )
  }

  ## Returning a ggplot object
  return(p)
}
