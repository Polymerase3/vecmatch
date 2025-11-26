#--divergent color palette vector-----------------------------------------------
vec_colors <- c(
  "#005b99", # Medium Dark Blue
  "#f1c232", # Light Golden Yellow
  "#d73a28", # Medium Dark Red
  "#00b25d", # Dark Green
  "#003d6e", # Dark Blue
  "#f57c00", # Darker Orange
  "#b51d14", # Dark Red
  "#1dd45c", # Bright Green
  "#042940", # Very Dark Blue
  "#f39c12", # Medium Orange
  "#e74c3c", # Bright Red
  "#009c49", # Medium Green
  "#00678a", # Deep Sky Blue
  "#a6c0f0", # Very Light Blue
  "#ddb310", # Golden Yellow
  "#6e8bce", # Light Blue
  "#ff7f00", # Bright Orange
  "#0084b4", # Bright Sky Blue
  "#d4af37", # Darker Gold
  "#4053d3" # Medium Blue
)

#--function to generate a list with colors from the vector----------------------
.generate_colors <- function(col_vector) {
  outlist <- list()
  i <- 1
  while (i <= length(col_vector)) {
    outlist <- append(outlist, list(col_vector[1:i]))
    i <- i + 1
  }
  return(outlist)
}

vec_col_list <- .generate_colors(vec_colors)
names(vec_col_list) <- c(
  "one", "two", "three", "four", "five",
  "six", "seven", "eight", "nine", "ten",
  "eleven", "twelve", "thirteen", "fourteen",
  "fifteen", "sixteen", "seventeen", "eighteen",
  "nineteen", "twenty"
)

#--function to generate palettes form the list of colors------------------------
# helper: parse "#RRGGBB" (or "#RGB") into numeric R,G,B in [0, 255]
.hex_to_rgb <- function(col) {
  col <- gsub("^#", "", col)

  # expand shorthand "#RGB" to "#RRGGBB"
  if (nchar(col) == 3L) {
    col <- paste0(
      substr(col, 1L, 1L), substr(col, 1L, 1L),
      substr(col, 2L, 2L), substr(col, 2L, 2L),
      substr(col, 3L, 3L), substr(col, 3L, 3L)
    )
  }

  c(
    R = strtoi(substr(col, 1L, 2L), 16L),
    G = strtoi(substr(col, 3L, 4L), 16L),
    B = strtoi(substr(col, 5L, 6L), 16L)
  )
}

# helper: clamp + convert back to "#RRGGBB"
.rgb_to_hex <- function(r, g, b) {
  r <- max(0, min(255, round(r)))
  g <- max(0, min(255, round(g)))
  b <- max(0, min(255, round(b)))
  sprintf("#%02X%02X%02X", r, g, b)
}

# simple linear interpolation in RGB space
.interpolate_palette <- function(cols, n) {
  m <- length(cols)

  if (m == 1L) {
    # only one anchor colour: just repeat it
    return(rep(cols, n))
  }

  # anchor positions in [0, 1]
  stops <- seq(0, 1, length.out = m)

  # matrix m x 3 of anchor RGB
  rgb_mat <- t(vapply(cols, .hex_to_rgb, numeric(3L)))

  # target positions
  pos <- if (n == 1L) 0.5 else seq(0, 1, length.out = n)

  out_rgb <- matrix(NA_real_, nrow = n, ncol = 3L)

  for (i in seq_len(n)) {
    x <- pos[i]

    if (x <= 0) {
      out_rgb[i, ] <- rgb_mat[1L, ]
    } else if (x >= 1) {
      out_rgb[i, ] <- rgb_mat[m, ]
    } else {
      # find interval [stops[j], stops[j+1]]
      j <- max(which(stops <= x))
      if (j >= m) j <- m - 1L

      t0 <- stops[j]
      t1 <- stops[j + 1L]
      alpha <- (x - t0) / (t1 - t0)

      out_rgb[i, ] <- (1 - alpha) * rgb_mat[j, ] + alpha * rgb_mat[j + 1L, ]
    }
  }

  vapply(
    seq_len(n),
    function(i) .rgb_to_hex(out_rgb[i, 1L], out_rgb[i, 2L], out_rgb[i, 3L]),
    character(1L)
  )
}

.generate_palettes <- function(n,
                               col_list = vec_col_list,
                               type = c("discrete", "continuous")) {
  type <- match.arg(type)
  cols <- col_list[[n]]

  palette <- switch(
    type,
    discrete   = cols,
    continuous = .interpolate_palette(cols, n)
  )

  structure(palette, name = names(col_list[[n]]), class = "palette")
}

#--ggplot2 functions to add scale-----------------------------------------------
scale_fill_vecmatch <- function(n, type = c("discrete", "continuous")) {
  if (type == "discrete") {
    ggplot2::scale_fill_manual(values = .generate_palettes(
      n = n,
      type = "discrete"
    ))
  } else {
    ggplot2::scale_fill_gradientn(colors = .generate_palettes(
      n = 3,
      type = "continuous"
    ))
  }
}

scale_color_vecmatch <- function(n, type = c("discrete", "continuous")) {
  if (type == "discrete") {
    ggplot2::scale_colour_manual(values = .generate_palettes(
      n = n,
      type = "discrete"
    ))
  } else {
    ggplot2::scale_colour_gradientn(colors = .generate_palettes(
      n = 3,
      type = "continuous"
    ))
  }
}

scale_colour_vecmatch <- scale_color_vecmatch

#--ggplot2 functions to add scale to mosaic plots
scale_x_product <- function(coords) {
  # subset data (only bottoms)
  coords_sub <- as.data.frame(subset(coords, coords$b == 0))

  # predefining vars
  label_position <- list()

  # calculating label positions and defining names
  for (i in seq_len(nrow(coords_sub))) {
    label_position$pos[i] <- (coords_sub$l[i] + coords_sub$r[i]) / 2
    label_position$name[i] <- as.character(coords_sub[i, 1])
  }

  # adding scale to plot
  ggplot2::scale_x_continuous(
    breaks = label_position$pos,
    labels = label_position$name
  )
}

scale_y_product <- function(coords) {
  # subset data (only lefts)
  coords_sub <- as.data.frame(subset(coords, coords$l == 0))

  if (!is.null(coords_sub$facet)) {
    coords_sub <- subset(
      coords_sub,
      coords_sub$facet == unique(coords_sub$facet)[1]
    )
    coords_sub <- coords_sub[, -which(colnames(coords_sub) == "facet")]
  }

  # Predefine vars
  label_position <- list()
  single_aes <- unique(coords_sub$level) == 1

  # Calculate positions and labels
  for (i in seq_len(nrow(coords_sub))) {
    sub_data <- coords_sub[i, ]

    # Ensure 'b' and 't' are treated as numeric
    pos <- (as.numeric(sub_data$b) + as.numeric(sub_data$t)) / 2
    label_position$pos[i] <- pos
    label_position$name[i] <- if (single_aes) {
      as.character(sub_data[, 1]) # First column for single aesthetics
    } else {
      as.character(sub_data[, which(names(sub_data) == "level") + 1])
    }
  }

  # adding scale to plot
  ggplot2::scale_y_continuous(
    breaks = label_position$pos,
    labels = label_position$name
  )
}

#--modified facet wrap to add scales to each plot
#-- solution by:
#-- https://dewey.dunnington.ca/post/2018/modifying-facet-scales-in-ggplot2/
#--simplified and modified
scale_facet <- function(which_facet, scale_arg) {
  structure(list(which_facet = which_facet, scale_arg = scale_arg),
    class = "scale_facet"
  )
}

FacetWrapScales <- ggplot2::ggproto("FacetWrapScales",
  ggplot2::FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    ## initial scales --> facet_wrap() deafult beahviour
    scales <- ggplot2::ggproto_parent(
      ggplot2::FacetWrap,
      self
    )$init_scales(
      layout,
      x_scale,
      y_scale,
      params
    )

    if (is.null(params$scales_custom)) {
      return(scales)
    }

    ## override the chosen scales based on scale_facet provided to facet_wrap
    for (scale_facet in params$scales_custom) {
      which_facet <- scale_facet$which_facet
      scale_arg <- scale_facet$scale_arg

      if ("x" %in% scale_arg$aesthetics) {
        if (!is.null(scales$x)) {
          scales$x[[which_facet]] <- scale_arg$clone()
        }
      } else if ("y" %in% scale_arg$aesthetics) {
        if (!is.null(scales$y)) {
          scales$y[[which_facet]] <- scale_arg$clone()
        }
      }
    }
    return(scales)
  }
)

facet_wrap_scales <- function(..., scales_custom = NULL) {
  facet_old <- ggplot2::facet_wrap(...)

  ## scales_custom handling
  if (inherits(scales_custom, "scale_facet")) {
    scales_custom <- list(scales_custom)
  }

  ## add scales_custom to params
  facet_old$params$scales_custom <- scales_custom

  ## create new facet_wrap
  ggplot2::ggproto(NULL, FacetWrapScales,
    shrink = facet_old$shrink,
    params = facet_old$params
  )
}

# define the replace function
"%+replace%" <- ggplot2::"%+replace%"

#' Fixing bug in productplots::prodcalc
#' @keywords internal
#' @return No return value, called for side effects
#' @export
vspine <- function(...) {
  productplots::vspine(...)
}

#' Fixing bug in productplots::prodcalc
#' @keywords internal
#' @return No return value, called for side effects
#' @export
hspine <- function(...) {
  productplots::hspine(...)
}

# Fixing hmisc dependency on the mean_cl function
mean_ci <- function(x) {
  mean <- mean(x)
  ci_lower <- mean - 1.96 * stats::sd(x) / sqrt(length(x))
  ci_upper <- mean + 1.96 * stats::sd(x) / sqrt(length(x))

  data.frame(y = mean, ymin = ci_lower, ymax = ci_upper)
}
