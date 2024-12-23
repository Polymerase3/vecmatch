## Defining new geom to raincloud plots based on:
## https://gist.github.com/benmarwick/2a1bb0133ff568cbe28d/

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show_legend = NA, inherit_aes = TRUE, ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show_legend,
    inherit.aes = inherit_aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}


GeomFlatViolin <-
  ggplot2::ggproto("GeomFlatViolin", ggplot2::Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (ggplot2::resolution(data$x, FALSE) * 0.9)

      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data <- split(data, data$group)

      data <- lapply(data, function(group_data) {
        group_data$ymin <- min(group_data$y)
        group_data$ymax <- max(group_data$y)
        group_data$xmin <- group_data$x
        group_data$xmax <- group_data$x + group_data$width / 2
        return(group_data)
      })

      data <- do.call(rbind, data)

      data
    },
    draw_group = function(data, panel_scales, coord) {
      # Find the points for the line to go all the way around
      data <- transform(data,
        xminv = x,
        xmaxv = x + violinwidth * (xmax - x)
      )

      # Make sure it's sorted properly to draw the outline
      # Create the first transformed dataset where x = xminv and arrange by y
      data1 <- transform(data, x = xminv)
      data1 <- data1[order(data1$y), ]

      # Create the second transformed dataset where x = xmaxv and arrange by -y
      # (descending order)
      data2 <- transform(data, x = xmaxv)
      data2 <- data2[order(-data2$y), ]

      # Combine both datasets using rbind
      newdata <- rbind(data1, data2)

      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1, ])

      ggplot2:::ggname(
        "geom_flat_violin",
        ggplot2::GeomPolygon$draw_panel(
          newdata,
          panel_scales,
          coord
        )
      )
    },
    draw_key = ggplot2::draw_key_polygon,
    default_aes = ggplot2::aes(
      weight = 1, colour = "grey20", fill = "white", linewidth = 0.5,
      alpha = NA, linetype = "solid"
    ),
    required_aes = c("x", "y")
  )
