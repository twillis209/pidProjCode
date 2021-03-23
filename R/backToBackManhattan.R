#' @title Writes out a back-to-back Manhattan plot.
#' 
#' @description Main back-to-back plot can be supplemented with a third plot and with gene structure diagrams.
#'
#' @details At present we require that the p-value column in a GenomicRanges column be labelled 'P' or 'p'.
#'
#' @param topGRanges GenomicRanges object containing data for top plot
#' @param bottomGRanges GenomicRanges object containing data for bottom plot
#' @param topLabel Axis label for top Manhattan plot of back-to-back pair
#' @param bottomLabel Axis label for bottom Manhattan plot of back-to-back pair
#' @param main Main title
#' @param ymax Maximum value for y axis
#' @param tickDist Distance in basepairs between chromosome diagram ticks
#' @param axisLabelMargin Size of margin between axis and its label
#' @param mainTitleCex Scaling factor for main title
#' @param axisLabelCex Scaling factor for axis label
#' @param plotGenes Flag to have karyoploteR plot gene structures beneath Manhattan
#' @param geneNamesCex Scaling factor for gene labels
#' @param chromNamesCex Scaling factor for chromosome labels
#' @param thirdGRanges GenomicRanges object containing data for third plot to be placed above back-to-back Manhattans
#' @param thirdLabel Axis label for third plot
#' @param chromosomes List of chromosomes to plot. Expects strings of the form 'chrx' where 'x' is from the set {1, 2, ..., 22, X, Y} 
#' @param zoom Coordinates for interval to magnify. Expects a string of the form 'chrx:a-b' where x is the chromosome number a and b are basepair coordinates.
#'
#' @importFrom karyoploteR plotKaryotype kpAddBaseNumbers kpAddChromosomeNames kpAddLabels kpAxis kpPlotGenes kpPlotManhattan makeGenesDataFromTxDb addGeneNames mergeTranscripts
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom GenomicRanges mcols mcols<-
#' @export
#' 
#' @examples
backToBackManhattan <- function(topGRanges, bottomGRanges, topLabel, bottomLabel, main, ymax=15, tickDist=1e5, axisLabelMargin=0.03, mainTitleCex=2.7, axisLabelCex=1.8, plotGenes=F, geneNamesCex=1, chromNamesCex=2, thirdGRanges=NULL, thirdLabel=NULL, chromosomes=NULL, zoom=NULL) {

  if(!is.null(chromosomes) & !is.null(zoom)) {
    stop("Cannot specify both \'chromosomes\' and \'zoom\'")
  }

  if(!is.null(thirdGRanges)) {
    if(is.null(thirdLabel)) {
      stop("Need to specify \'thirdLabel\' if specifying \'thirdGRanges\'")
      }
    if(!('P' %in% names(mcols(thirdGRanges)) | 'p' %in% names(mcols(thirdGRanges)))) {
      stop('Need a p-value column labelled \'P\' or \'p\' in thirdGRanges')
    } 
  }

  if(plotGenes & is.null(zoom)) {
    stop("Can only plot gene tracks if \'zoom\' coordinates are specified")
  }

  if(!('P' %in% names(mcols(topGRanges)) | 'p' %in% names(mcols(topGRanges)))) {
    stop('Need a p-value column labelled \'P\' or \'p\' in topGRanges')
  }

  if(!('P' %in% names(mcols(bottomGRanges)) | 'p' %in% names(mcols(bottomGRanges)))) {
    stop('Need a p-value column labelled \'P\' or \'p\' in bottomGRanges')
  } 

  if(!is.null(chromosomes)) {
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL, chromosomes=chromosomes)
    kpAddBaseNumbers(kp, add.units=T, cex=1, tick.dist=tickDist)
  } else if(!is.null(zoom)) {
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL, zoom=zoom)
    kpAddBaseNumbers(kp, add.units=T, cex=1, tick.dist=tickDist)
  } else {
    kp <- plotKaryotype(plot.type=4, labels.plotter=NULL)
    kpAddChromosomeNames(kp, col='black',srt=90,cex=chromNamesCex)
  } 

  title(main=main, cex.main= mainTitleCex)

  if(is.null(thirdGRanges)) {
    if(plotGenes) {
      # two tracks, plot genes
      kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.6, r1=1, cex=axisLabelCex, label.margin = axisLabelMargin)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.6)
      kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.6, r1=1, ymax=ymax)

      kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0.6, r1=0.2, cex=axisLabelCex, label.margin = axisLabelMargin)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.6, r1=0.2)
      kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.6, r1=0.2, ymax=ymax, points.col = "2blues")
    } else {
      # two tracks, do not plot genes
      kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.5, r1=1, cex=axisLabelCex, label.margin = axisLabelMargin)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.5)
      kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.5, r1=1, ymax=ymax)

      kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0, r1=0.5, cex=axisLabelCex, label.margin = axisLabelMargin)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.5, r1=0)
      kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.5, r1=0, ymax=ymax, points.col = "2blues")
    }
  } else {
    if(plotGenes) {
      # three-track Manhattan, plot genes
      kpAddLabels(kp, labels = thirdLabel, srt=90, pos=3, r0=0.74, r1=1, cex=axisLabelCex, label.margin = axisLabelMargin)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.74, r1=1)
      kp <- kpPlotManhattan(kp, data=thirdGRanges, r0=0.74, r1=1, ymax=ymax)

      kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.46, r1=0.72, cex=axisLabelCex, label.margin = axisLabelMargin)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.46, r1=0.72)
      kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.46, r1=0.72, ymax=ymax)

      kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0.2, r1=0.46, cex=axisLabelCex, label.margin = axisLabelMargin)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.46, r1=0.2)
      kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.46, r1=0.2, ymax=ymax, points.col = "2blues")
    } else {
      # Three-track Manhattan
      kpAddLabels(kp, labels = thirdLabel, srt=90, pos=3, r0=0.68, r1=1, cex=axisLabelCex, label.margin = axisLabelMargin)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.68, r1=1)
      kp <- kpPlotManhattan(kp, data=thirdGRanges, r0=0.68, r1=1, ymax=ymax)

      kpAddLabels(kp, labels = topLabel, srt=90, pos=3, r0=0.33, r1=0.66, cex=axisLabelCex, label.margin = axisLabelMargin)
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.33, r1=0.66)
      kp <- kpPlotManhattan(kp, data=topGRanges, r0=0.33, r1=0.66, ymax=ymax)

      kpAddLabels(kp, labels = bottomLabel, srt=90, pos=3, r0=0, r1=0.33, cex=axisLabelCex, label.margin = axisLabelMargin)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kpAxis(kp, ymin=0, ymax=ymax, r0=0.33, r1=0)
      kp <- kpPlotManhattan(kp, data=bottomGRanges, r0=0.33, r1=0, ymax=ymax, points.col = "2blues")
      }
  }

  if(plotGenes) {
    genes.data<-makeGenesDataFromTxDb(karyoplot=kp, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)

    genes.data<-addGeneNames(genes.data)
    genes.data.merged<-mergeTranscripts(genes.data)
    kp<-kpPlotGenes(kp, data=genes.data.merged, r0=0, r1=0.2, cex=geneNameCex, gene.name.position='left')
  }
}
