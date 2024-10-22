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
.generate_palettes <- function(n, col_list = vec_col_list,
                               type = c("discrete", "continuous")) {
  cols <- col_list[[n]]
  type <- match.arg(type)

  palette <- switch(type,
    discrete = cols,
    continuous = grDevices::colorRampPalette(cols)(n)
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
