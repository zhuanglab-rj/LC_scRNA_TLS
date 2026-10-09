#!/usr/bin/env Rscript
# =============================================================================
# Custom ggplot2 geom for split violin plots
# Used in Figure 1E for IHC density comparison
#
# Original implementation by Jan Gleixner (@jan-glx)
# Adjustments by Wouter van der Bijl (@Axeman)
# =============================================================================

#' GeomSplitViolin - A ggproto object for split violin plots
#'
#' Draws violin plots split by group, allowing comparison of distributions
#' side by side within the same x-axis position.
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data,
      xminv = x - violinwidth * (x - xmin),
      xmaxv = x + violinwidth * (xmax - x)
    )

    grp <- data[1, "group"]
    newdata <- plyr::arrange(
      transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
      if (grp %% 2 == 1) y else -y
    )
    newdata <- rbind(
      newdata[1, ], newdata,
      newdata[nrow(newdata), ], newdata[1, ]
    )
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- create_quantile_segment_frame(
        data, draw_quantiles,
        split = TRUE, grp = grp
      )
      aesthetics <- data[
        rep(1, nrow(quantiles)),
        setdiff(names(data), c("x", "y")),
        drop = FALSE
      ]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
        grid::grobTree(
          GeomPolygon$draw_panel(newdata, ...),
          quantile_grob
        )
      )
    } else {
      ggplot2:::ggname("geom_split_violin",
        GeomPolygon$draw_panel(newdata, ...)
      )
    }
  }
)

#' Create quantile segment frame for split violin plots
#'
#' @param data Data frame with violin plot data
#' @param draw_quantiles Numeric vector of quantiles to draw
#' @param split Logical, whether to split the violin
#' @param grp Group identifier
#' @return Data frame with quantile segment coordinates
create_quantile_segment_frame <- function(data, draw_quantiles,
                                          split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)

  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2),
      group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2),
      group = rep(ys, each = 2)
    )
  }
}

#' geom_split_violin - Create split violin plots
#'
#' @param mapping Aesthetic mapping
#' @param data Data frame
#' @param stat Statistical transformation
#' @param position Position adjustment
#' @param ... Additional parameters
#' @param draw_quantiles Quantiles to draw
#' @param trim Whether to trim the violin
#' @param scale Scaling method
#' @param na.rm Remove NA values
#' @param show.legend Show legend
#' @param inherit.aes Inherit aesthetics
#' @return A ggplot2 layer
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., draw_quantiles = NULL,
                              trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      trim = trim, scale = scale, draw_quantiles = draw_quantiles,
      na.rm = na.rm, ...
    )
  )
}
