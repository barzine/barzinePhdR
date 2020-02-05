#' Theme Legacy theme grey
#'
#' Theme resetting the plot title at the center (as it was in the earlier version of ggplot2)
#' @import ggplot2
#' @importFrom ggplot2 element_text element_blank element_line
#' @inheritParams ggplot2::theme_grey
#' @return An object of class \code{\link[ggplot2]{theme}()}.
#' @export
#' @family themes legacy
#' @examples ggplot(ggplot2::mpg, aes(cyl, hwy)) + geom_point() + ggtitle('Centered title') + theme_legacy()
theme_legacy<-function (base_size = 11, base_family = "", base_line_size = base_size/22,
                        base_rect_size = base_size/22) {
  ret<-ggplot2::theme_grey()+ggplot2::theme(plot.title=element_text(hjust=0.5))
  return(ret)
}


#' Theme minimaliste: customised theme for ggplot2
#'
#' @param base_size Default: 11
#' @param base_family Default ""
#' @param ticks Boolean; whether the ticks should be drawn on the axes
#'
#' @inheritParams ggplot2::theme_bw
#' @return An object of class \code{\link[ggplot2]{theme}()}.
#' @export
#' @family themes minimaliste
#' @examples ggplot(ggplot2::mpg, aes(cyl, hwy)) + geom_point() + ggtitle('Centered title') + theme_minimaliste()
theme_minimaliste<-function(base_size = 11 , base_family = "", ticks = TRUE){
    ret <- ggplot2::theme_bw(base_family = base_family, base_size = base_size) +
      ggplot2::theme(panel.background = element_blank(),# panel.grid = element_blank(),
              panel.grid.major = element_line(colour = "gray30",size=0.02),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), #strip.background = element_blank(),
              plot.background = element_blank(),
              legend.background = element_blank(), legend.key = element_blank(),
              axis.line = element_blank(),
              plot.title=element_text(hjust=0.5)
              )
    if (!ticks) {
        ret <- ret + ggplot2::theme(axis.ticks = element_blank())
    }
    return(ret)
}

#' Theme black and white with centered title
#'
#' Theme resetting the plot title at the center (as it was in the earlier version of ggplot2)
#'
#' @inheritParams ggplot2::theme_bw
#' @return An object of class \code{\link[ggplot2]{theme}()}.
#' @export
#' @family themes cbw
#' @examples ggplot(ggplot2::mpg, aes(cyl, hwy)) + geom_point() + ggtitle('Centered title') + theme_cbw()
theme_cbw<-function(base_size = 11, base_family = "", base_line_size = base_size/22,
                             base_rect_size = base_size/22){
  ret <- ggplot2::theme_bw()+ggplot2::theme(plot.title=element_text(hjust=0.5))
  return(ret)
}

#' Extract the legend from a ggplot
#'
#' @description function developed by Hadley Wickham
#' found in this (by Michael Kuhn) blog:http://blog.mckuhn.de/2013/04/2d-plot-with-histograms-for-each.html
#' and was referencing: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#' @author Hadley Wickham
#'
#' @param a_ggplot a ggplot
#'
#' @return the legend
#' @export
#'
g_legend<-function(a_ggplot){
  tmp <- ggplot_gtable(ggplot_build(a_ggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' Create a palette with ggplot2 default colors
#'
#' @author John Colby \url{https://stackoverflow.com/a/8197703/2116422}
#'
#' @param n integer. Number of colours in the palette.
#'
#' @return vector of character strings corresponding to the palette colour hex codes.
#' @export
#'
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
