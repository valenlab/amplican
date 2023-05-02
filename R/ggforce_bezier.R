#' Create quadratic or cubic bezier curves [copied from ggforce]
#'
#' This set of functionality is copied from ggforce package due to dependency
#' issues on Bioconductor and is used internally (not exported) only.
#' This set of geoms makes it possible to connect points creating either
#' quadratic or cubic beziers. bezier works by calculating
#' points along the bezier and connecting these to draw the curve.
#'
#' @details
#' Input data is understood as a sequence of data points the first being the
#' start point, then followed by one or two control points and then the end
#' point. More than 4 and less than 3 points per group will throw an error.
#'
#' @section Aesthetics:
#' geom_link, geom_link2 and geom_lin0 understand the following aesthetics
#' (required aesthetics are in bold):
#'
#' - **x**
#' - **y**
#' - color
#' - size
#' - linetype
#' - alpha
#' - lineend
#'
#'
#' @section Computed variables:
#'
#' \describe{
#'  \item{x, y}{The interpolated point coordinates}
#'  \item{index}{The progression along the interpolation mapped between 0 and 1}
#' }
#'
#' @inheritParams ggplot2::geom_path
#' @inheritParams ggplot2::stat_identity
#'
#' @param n The number of points to create for each segment
#'
#' @author Thomas Lin Pedersen
#'
#' @name geom_bezier
#' @rdname geom_bezier
#'
#' @examples
#' beziers <- data.frame(
#'     x = c(1, 2, 3, 4, 4, 6, 6),
#'     y = c(0, 2, 0, 0, 2, 2, 0),
#'     type = rep(c('cubic', 'quadratic'), c(3, 4)),
#'     point = c('end', 'control', 'end', 'end', 'control', 'control', 'end')
#' )
#' help_lines <- data.frame(
#'     x = c(1, 3, 4, 6),
#'     xend = c(2, 2, 4, 6),
#'     y = 0,
#'     yend = 2
#' )
#' ggplot2::ggplot() + ggplot2::geom_segment(
#'   ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
#'   data = help_lines,
#'   arrow = ggplot2::arrow(length = ggplot2::unit(c(0, 0, 0.5, 0.5), 'cm')),
#'   colour = 'grey') +
#'   amplican:::geom_bezier(ggplot2::aes(x= x, y = y, group = type, linetype = type),
#'               data = beziers) +
#'   ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = point), data = beziers)
#'
NULL

#' @importFrom ggplot2 ggproto Stat
#'
StatBezier <- ggproto(
  'StatBezier', Stat, compute_layer = function(self, data, params, panels) {
    if (is.null(data)) return(data)
    nControls <- table(data$group)
    controlRange <- range(nControls)
    if (min(controlRange) < 3 || max(controlRange) > 4) {
      stop('Only support for quadratic and cubic beziers')
    }
    data <- data[order(data$group),]
    paths <- getBeziers(data$x, data$y, data$group, params$n)
    paths <- data.frame(x = paths$paths[,1], y = paths$paths[,2],
                        group = paths$pathID)
    paths$index <- rep(seq(0, 1, length.out = params$n), length(nControls))
    dataIndex <- rep(match(unique(data$group), data$group), each = params$n)
    cbind(paths, data[dataIndex, !names(data) %in% c('x', 'y', 'group'),
                      drop = FALSE])
  },
  required_aes = c('x', 'y'),
  extra_params = c('na.rm', 'n')
)

#' @rdname geom_bezier
#' @importFrom ggplot2 layer
#'
stat_bezier <- function(mapping = NULL, data = NULL, geom = "path",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        n = 100, inherit.aes = TRUE, ...) {
  layer(
    stat = StatBezier, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, n = n, ...)
  )
}

#' @rdname geom_bezier
#' @importFrom ggplot2 layer
#'
geom_bezier <- function(mapping = NULL, data = NULL, stat = "bezier",
                        position = "identity", arrow = NULL, lineend = "butt",
                        na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                        n = 100, ...) {
  layer(data = data, mapping = mapping, stat = stat,
        geom = ggplot2::GeomPath,
        position = position, show.legend = show.legend,
        inherit.aes = inherit.aes, params = list(
          arrow = arrow, lineend = lineend, na.rm = na.rm, n = n, ...))
}
