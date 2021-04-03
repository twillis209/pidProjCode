#' @title Writes out a stratified Q-Q plot.
#'
#' @description
#'
#' @details Note that this function does not do the heavy lifting of styling the plot's aesthetics.
#' 
#' @param dataFrame \code{data.frame} containing p-values and conditional p-values
#' @param principalValueLabel
#' @param conditionalValueLabel
#' @param thresholds
#'
#' @importFrom ggplot2 ggplot geom_line geom_abline
#' @export
#' TODO write unconditional Q-Q plot example
#' TODO write conditional Q-Q plot example, include modification
#' TODO handle marginal Q-Q plot
#' @examples
#'
#' sample_data <- data.frame(p=c(runif(9e5), pnorm(rnorm(1e5))), q=c(runif(9e5), pnorm(rnorm(1e5))))
#'
#' stratified_qqplot(sample_data, 'p', 'q')
#' 
stratified_qqplot <- function(dataFrame, principalValueLabel, conditionalValueLabel = NULL, thresholds = c(1, 1e-1, 1e-2, 1e-3, 1e-4)) {

  dataFrame$negLogP <- -log10(dataFrame[, c(principalValueLabel)])

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
