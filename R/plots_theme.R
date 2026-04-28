#' ProteoForge ggplot2 theme
#'
#' A minimalist white-background theme using the Okabe-Ito colorblind-safe
#' palette. Suitable for publication-quality plots and grayscale printing.
#'
#' @param base_size   Base font size. Default 11.
#' @param base_family Base font family. Default `""` (device default).
#'
#' @return A `ggplot2::theme` object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(hp, mpg)) + geom_point() + theme_proteoforge()
#'
#' @export
theme_proteoforge <- function(base_size = 11, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major  = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(colour = "grey60", fill = NA,
                                                linewidth = 0.6),
      axis.line         = ggplot2::element_blank(),
      axis.ticks        = ggplot2::element_line(colour = "grey60", linewidth = 0.4),
      strip.background  = ggplot2::element_rect(fill = "grey95", colour = "grey60",
                                                linewidth = 0.4),
      strip.text        = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.key        = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(size = ggplot2::rel(1.1),
                                                face = "bold"),
      plot.subtitle     = ggplot2::element_text(colour = "grey40",
                                                size = ggplot2::rel(0.9)),
      plot.caption      = ggplot2::element_text(colour = "grey50",
                                                size = ggplot2::rel(0.75)),
      plot.background   = ggplot2::element_blank(),
      complete          = FALSE
    )
}

#' Okabe-Ito colorblind-safe palette
#'
#' Returns n colours from the Okabe-Ito palette, which is distinguishable for
#' the most common forms of colour-blindness and in grayscale print.
#'
#' @param n Number of colours to return. Must be ≤ 8.
#'
#' @return Character vector of hex colour codes.
#'
#' @examples
#' pf_palette(4)
#'
#' @export
pf_palette <- function(n = 8L) {
  oi <- c("#E69F00","#56B4E9","#009E73","#F0E442",
          "#0072B2","#D55E00","#CC79A7","#000000")
  if (n > length(oi)) {
    warning("Requested ", n, " colours but palette only has ", length(oi),
            ". Colours will be recycled.", call. = FALSE)
    oi <- rep(oi, ceiling(n / length(oi)))
  }
  oi[seq_len(n)]
}

#' Continuous fill scale using the ProteoForge diverging palette
#'
#' @param ... Additional arguments passed to `ggplot2::scale_fill_gradientn`.
#'
#' @return A ggplot2 scale object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg, fill = hp)) +
#'   geom_point(shape = 21, size = 3) + scale_fill_pf()
#'
#' @export
scale_fill_pf <- function(...) {
  ggplot2::scale_fill_gradientn(
    colours = c("#0072B2","white","#D55E00"),
    ...
  )
}

#' Discrete colour scale using the Okabe-Ito palette
#'
#' @param ... Additional arguments passed to `ggplot2::scale_colour_manual`.
#'
#' @return A ggplot2 scale object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(hp, mpg, colour = factor(cyl))) +
#'   geom_point() + scale_colour_pf()
#'
#' @export
scale_colour_pf <- function(...) {
  ggplot2::scale_colour_manual(values = pf_palette(8), ...)
}

#' @rdname scale_colour_pf
#' @export
scale_color_pf <- scale_colour_pf
