#' @title ggplot2 theme for Zoom presentation
#'
#' @import ggplot2
#' @export
zoom_theme <- function() {
  theme_bw()+theme(
              axis.title = element_text(size=15),
              plot.title=element_text(hjust=0.5, size=15),
              strip.text=element_text(size=8),
              axis.text.x=element_text(size=8, angle=30, color="black"),
              axis.text.y=element_text(size=12, color="black"),
              legend.title=element_text(size=12),
              legend.text=element_text(size=12)
            )
}

#' @title ggplot2 theme for first-year report
#'
#' @import ggplot2
#' @export
fyr_theme <- function() {
  theme_bw()+theme(
               axis.title = element_text(size=12),
               plot.title=element_text(hjust=0.5, size=12),
               strip.text=element_text(size=10),
               axis.text.x=element_text(size=10, angle=30, color="black"),
               axis.text.y=element_text(size=10, color="black"),
               legend.title=element_text(size=10),
               legend.text=element_text(size=10)
             )
}
