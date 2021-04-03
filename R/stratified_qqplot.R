#' @title Writes out a stratified Q-Q plot.
#'
#' @description
#'
#' @details Note that this function does not do the heavy lifting of styling the plot's aesthetics.
#' 
#' @param dataFrame \code{data.frame} containing p-values and conditional p-values
#' @param principalValueLabel label of principal p-value column in \code{dataFrame}
#' @param conditionalValueLabel label of conditional trait column in \code{dataFrame}
#' @param thresholds threshold values to define strata
#'
#' @import ggplot2 
#'
#' @return ggplot object 
#' @export
#' 
stratified_qqplot <- function(dataFrame, principalValueLabel, conditionalValueLabel = NULL, thresholds = c(1, 1e-1, 1e-2, 1e-3, 1e-4)) {

  dataFrame$negLogP <- -log10(dataFrame[, principalValueLabel])

  if(is.null(conditionalValueLabel)) {
    daf <- dataFrame[, c(principalValueLabel, 'negLogP')]
    daf <- daf[order(daf[,principalValueLabel]), ]
    daf$pp <- -log10(ppoints(nrow(daf)))
    daf$threshold <- factor(c(1))
  } else {
    dataFrame <- dataFrame[, c(principalValueLabel, conditionalValueLabel, 'negLogP')]

    dafs <- list()

    for(i in seq_along(thresholds)) {
      daf <- subset(dataFrame, get(conditionalValueLabel) < thresholds[i])
      daf <- daf[order(daf[ , principalValueLabel]) , ]
      daf$pp <- -log10(ppoints(nrow(daf)))
      daf$threshold <- factor(thresholds)[i]
      dafs[[i]] <- daf
    }

      daf <- do.call(rbind, dafs)
    }

  ggplot(data=daf) + geom_line(aes(x = pp, y = negLogP, group = threshold, colour = threshold)) +geom_abline(intercept=0,slope=1, linetype="dashed")
}
