#' Examine the Distribution of Continuous Data with Raincloud Plots
#'
#' @description
#' The `raincloud()` function allows to generate raincloud plots for continuous
#' data in an easy and uncomplicated way. Raincloud plots consist of three main
#' elements:
#' /item{Distribution plots}{For example density plots, histograms or violin
#' plots wit the mean values of respective groups}
#' /item{Jittered point plots}{Depicting the underlying distribution of the data
#' in the rawest form}
#' /item{Boxplots}
#' The function is based on the `ggplot2` package, which must already be
#' preinstalled
#'
#' @param data A non-empty object of the class `data.frame` with at least one
#' numeric column
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
                      significance = FALSE,        ## not functional
                      limits = NULL,
                      jitter = 0.1,
                      alpha = 0.4,
                      save = FALSE,
                      plot.name = NULL,
                      overwrite = FALSE) {        ## not functional
  ############################ INPUT CHECKING###################################
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

  ## Check character and valid name for plot.name
  if(!is.null(plot.name)) {
    chk::chk_character(plot.name)
    .check_extension(plot.name, x_name = 'plot.name',
                     ext_vec = c('.png', '.PNG', '.pdf', '.PDF'))
  }

  if(save == TRUE && is.null(plot.name)) {
    chk::abort_chk('If save == TRUE then the `plot.name` argument has to be specified.')
  }

  ## Check logical for overwrite
  chk::chk_logical(overwrite)

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
  #--defining necessary variables-----------------------------------------------
  "%+replace%" <- ggplot2::"%+replace%"
  rain_height <- 0.1

  ## Unique values in grouping variables (necessary to define the palette)
  pal_len <- length(unique(data[, symlist[["group"]]]))
  if (is.null(symlist[["group"]]) || pal_len == 0) pal_len <- 1

  ##--defining the main ggplot formula------------------------------------------
  main <- paste0(
    "ggplot2::ggplot(data, ggplot2::aes(x = '', y = y",
    ifelse(is.null(symlist[["group"]]), "))", paste0(
      ", fill = ",
      symlist[["group"]], ", color = ", symlist[["group"]],
      "))"
    ))
  )

  ##--defining the geom_jitter
  geom_jitter <- paste0(
    'ggplot2::geom_jitter( size = 2, alpha = alpha, show.legend = FALSE, ',
    ifelse(pal_len == 1,
           'position = ggplot2::position_jitter(width = jitter))',
           'position = ggplot2::position_jitterdodge(jitter.width = jitter,
           dodge.width = 0.25))')
  )

  ##--defining the geom_boxplot
  geom_boxplot <- paste0(
    'ggplot2::geom_boxplot(width = 0.1, alpha = alpha, show.legend = FALSE, ',
    ifelse(pal_len == 1,
           'position = ggplot2::position_nudge(x = -0.22))',
           'position = ggpp::position_dodgenudge(width = 0.2, x = -0.22))')
  )

  ##--defining the stat_summary
  stat_summ <- paste0(
    'ggplot2::stat_summary(fun.data = ggplot2::mean_cl_normal,
                           show.legend = FALSE, ',
    ifelse(pal_len == 1,
           'position = ggplot2::position_nudge(x = rain_height * 3))',
           'position = ggpp::position_dodgenudge(x = rain_height * 3, width = 0.1))')
  )

  #--defining the ggplot object-------------------------------------------------
  p <- eval(parse(text = main)) +
    ## halfs of the violin plots
    geom_flat_violin(
      trim = FALSE, alpha = alpha,
      position = ggplot2::position_nudge(x = rain_height + 0.05)
    ) +
    ## datapoints
    eval(parse(text = geom_jitter)) +
    ## boxplots
    eval(parse(text = geom_boxplot)) +
    ## stat_summary
    eval(parse(text = stat_summ)) +
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
    suppressMessages(ggplot2::ggsave(plot.name,
           plot = p, dpi = 300, create.dir = TRUE)
    )
  }

  ## Returning a ggplot object
  return(p)
}
