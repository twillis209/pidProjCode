#' @title Writes out a back-to-back Manhattan plot.
#' 
#' @description 
#'
#' @details 
#'
#' @param topGRanges GenomicRanges object containing data for top plot
#' @param bottomGRanges GenomicRanges object containing data for bottom plot
#' @param outputFile
#' @param topLabel
#' @param bottomLabel
#' @param main
#' @param width
#' @param height
#' @param ymax
#' @param thirdGRanges
#' @param thirdLabel
#' @param chromosomes
#' @param zoom
#'
#' @importFrom karyoploteR plotKaryotype kpAddBaseNumbers kpAddChromosomeNames kpAddLabels kpAxis kpPlotGenes kpPlotManhattan
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @export
#'
#' @examples
#' 
#' 
backToBackManhattan <- function(topGRanges, bottomGRanges, outputFile, topLabel, bottomLabel, main, width=1800, height=1000, ymax=15, thirdGRanges=NULL, thirdLabel=NULL, chromosomes=NULL, zoom=NULL) {

  if(!is.null(chromosomes) & !is.null(zoom)) {
    stop("Cannot specify both \'chromosomes\' and \'zoom\'")
  }

  if(!is.null(thirdGRanges)) {
    if(is.null(thirdLabel)) stop("Need to specify \'thirdLabel\' if specifying \'thirdGRanges\'")
  }

  png(outputFile, width=width, height=height)

  if(!is.null(chromosomes)) {
    # TODO compute sensible tick.dist
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL, chromosomes=chromosomes)
    kpAddBaseNumbers(kp, add.units=T, cex=1, tick.dist=5e6)
  } else if(!is.null(zoom)) {
    # TODO compute sensible tick.dist
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL, zoom=zoom)
    kpAddBaseNumbers(kp, add.units=T, cex=1, tick.dist=1e5)
  } else {
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL)
    kpAddChromosomeNames(kp, col='black',srt=90,cex=2)
  } 

  title(main=main, cex.main=2.7)

  if(is.null(thirdGRanges)) {
    kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.5, r1=1, cex=1.8, label.margin = 0.025)
    kpAxis(kp, ymin=0, ymax=ymax, r0=0.5)
    kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.5, r1=1, ymax=ymax)

    kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0, r1=0.5, cex=1.8, label.margin = 0.025)
    kpAxis(kp, ymin=0, ymax=ymax, r0=0.5, r1=0)
    kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.5, r1=0, ymax=ymax, points.col = "2blues")
  } else {
    # Three-track Manhattan
    kpAddLabels(kp, labels = thirdLabel, srt=90, pos=3, r0=0.68, r1=1, cex=1.8, label.margin = 0.025)
    kpAxis(kp, ymin=0, ymax=ymax, r0=0.68, r1=1)
    kp <- kpPlotManhattan(kp, data=pidGRanges, r0=0.68, r1=1, ymax=ymax)

    kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.33, r1=0.66, cex=1.8, label.margin = 0.025)
    kpAxis(kp, ymin=0, ymax=ymax, r0=0.33, r1=0.66)
    kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.33, r1=0.66, ymax=ymax)

    kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0, r1=0.33, cex=1.8, label.margin = 0.025)
    # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
    kpAxis(kp, ymin=0, ymax=ymax, r0=0.33, r1=0)
    kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.33, r1=0, ymax=ymax, points.col = "2blues")
  }

  dev.off()
}
