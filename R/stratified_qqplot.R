#' @title Writes out a stratified Q-Q plot.
#'
#' @description
#'
#' @details Note that this function does not do the heavy lifting of styling the plot's aesthetics.
#' 
#' @param data_frame \code{data.frame} containing p-values and conditional p-values
#' @param prin_value_label label of principal p-value column in \code{data_frame}
#' @param cond_value_label label of conditional trait column in \code{data_frame}
#' @param thresholds threshold values to define strata
#'
#' @import ggplot2 
#'
#' @return ggplot object 
#' @export
#' 
stratified_qqplot <- function(data_frame, prin_value_label, cond_value_label = NULL, thresholds = c(1, 1e-1, 1e-2, 1e-3, 1e-4)) {

  data_frame$negLogP <- -log10(data_frame[, prin_value_label])

  if(is.null(cond_value_label)) {
    daf <- data_frame[, c(prin_value_label, 'negLogP')]
    daf <- daf[order(daf[,prin_value_label]), ]
    daf$pp <- -log10(ppoints(nrow(daf)))
    daf$threshold <- factor(c(1))
  } else {
    data_frame <- data_frame[, c(prin_value_label, cond_value_label, 'negLogP')]

    dafs <- list()

    for(i in seq_along(thresholds)) {
      daf <- subset(data_frame, get(cond_value_label) < thresholds[i])
      daf <- daf[order(daf[ , prin_value_label]) , ]
      daf$pp <- -log10(ppoints(nrow(daf)))
      daf$threshold <- factor(thresholds)[i]
      dafs[[i]] <- daf
    }

      daf <- do.call(rbind, dafs)
    }

  ggplot(data=daf) + geom_line(aes(x = pp, y = negLogP, group = threshold, colour = threshold)) +geom_abline(intercept=0,slope=1, linetype="dashed")
}
