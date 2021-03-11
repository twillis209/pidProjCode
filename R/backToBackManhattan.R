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
#' @param tickDist
#' @param plotGenes  
#' @param thirdGRanges
#' @param thirdLabel
#' @param chromosomes
#' @param zoom
#'
#' @importFrom karyoploteR plotKaryotype kpAddBaseNumbers kpAddChromosomeNames kpAddLabels kpAxis kpPlotGenes kpPlotManhattan makeGenesDataFromTxDb addGeneNames mergeTranscripts
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @export
#' Use cases: - two plots: c(0, 0.5, 1.0)
#'            - two plots and a gene track: (0, 0.2), (0.2, 0.6), (0.6, 1.0)
#'            - three plots: (0, 0.33), (0.33, 0.66), (0.68, 1.0)
#'            - three plots and a gene track: (0, 0.2), (0.2, 0.46), (0.46, 0.72), (0.74, 1.0)
#' @examples
backToBackManhattan <- function(topGRanges, bottomGRanges, topLabel, bottomLabel, main, outputFile=NULL, width=1800, height=1000, ymax=15, tickDist=1e5, plotGenes=F, thirdGRanges=NULL, thirdLabel=NULL, chromosomes=NULL, zoom=NULL) {

  if(!is.null(chromosomes) & !is.null(zoom)) {
    stop("Cannot specify both \'chromosomes\' and \'zoom\'")
  }

  if(!is.null(thirdGRanges)) {
    if(is.null(thirdLabel)) stop("Need to specify \'thirdLabel\' if specifying \'thirdGRanges\'")
  }

  if(plotGenes & is.null(zoom)) {
    stop("Can only plot gene tracks if \'zoom\' coordinates are specified")
  }

  if(!is.null(outputFile)) {
    png(outputFile, width=width, height=height)
  }

  if(!is.null(chromosomes)) {
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL, chromosomes=chromosomes)
    kpAddBaseNumbers(kp, add.units=T, cex=1, tick.dist=tickDist)
  } else if(!is.null(zoom)) {
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL, zoom=zoom)
    kpAddBaseNumbers(kp, add.units=T, cex=1, tick.dist=tickDist)
  } else {
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL)
    kpAddChromosomeNames(kp, col='black',srt=90,cex=2)
  } 

  title(main=main, cex.main=2.7)

  if(is.null(thirdGRanges)) {
    if(plotGenes) {
      # two tracks, plot genes
      kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.6, r1=1, cex=1.8, label.margin = 0.025)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.6)
      kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.6, r1=1, ymax=ymax)

      kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0.6, r1=0.2, cex=1.8, label.margin = 0.025)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.6, r1=0.2)
      kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.6, r1=0.2, ymax=ymax, points.col = "2blues")
    } else {
      # two tracks, do not plot genes
      kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.5, r1=1, cex=1.8, label.margin = 0.025)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.5)
      kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.5, r1=1, ymax=ymax)

      kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0, r1=0.5, cex=1.8, label.margin = 0.025)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.5, r1=0)
      kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.5, r1=0, ymax=ymax, points.col = "2blues")
    }
  } else {
    if(plotGenes) {
      # three-track Manhattan, plot genes
      kpAddLabels(kp, labels = thirdLabel, srt=90, pos=3, r0=0.74, r1=1, cex=1.8, label.margin = 0.025)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.74, r1=1)
      kp <- kpPlotManhattan(kp, data=thirdGRanges, r0=0.74, r1=1, ymax=ymax)

      kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.46, r1=0.72, cex=1.8, label.margin = 0.025)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.46, r1=0.72)
      kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.46, r1=0.72, ymax=ymax)

      kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0.2, r1=0.46, cex=1.8, label.margin = 0.025)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.46, r1=0.2)
      kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.46, r1=0.2, ymax=ymax, points.col = "2blues")
    } else {
      # Three-track Manhattan
      kpAddLabels(kp, labels = thirdLabel, srt=90, pos=3, r0=0.68, r1=1, cex=1.8, label.margin = 0.025)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.68, r1=1)
      kp <- kpPlotManhattan(kp, data=thirdGRanges, r0=0.68, r1=1, ymax=ymax)

      kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.33, r1=0.66, cex=1.8, label.margin = 0.025)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.33, r1=0.66)
      kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.33, r1=0.66, ymax=ymax)

      kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0, r1=0.33, cex=1.8, label.margin = 0.025)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.33, r1=0)
      kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.33, r1=0, ymax=ymax, points.col = "2blues")
      }
  }

  if(plotGenes) {
    genes.data<-makeGenesDataFromTxDb(karyoplot=kp, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene)

    genes.data<-addGeneNames(genes.data)
    genes.data.merged<-mergeTranscripts(genes.data)
    kp<-kpPlotGenes(kp, data=genes.data.merged, r0=0, r1=0.2, cex=1.0, gene.name.position='left')
  }

  if(!is.null(outputFile)) {
    dev.off()
  }
}
