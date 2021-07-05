#' @title Writes out a back-to-back Manhattan plot.
#'
#' @description Main back-to-back plot can be supplemented with a third plot and with gene structure diagrams.
#'
#' @details At present we require that the p-value column in a GenomicRanges column be labelled 'P' or 'p'.
#'
#' @param top_gRanges GenomicRanges object containing data for top plot
#' @param bottom_gRanges GenomicRanges object containing data for bottom plot
#' @param top_label Axis label for top Manhattan plot of back-to-back pair
#' @param bottom_label Axis label for bottom Manhattan plot of back-to-back pair
#' @param main Main title
#' @param ymax Maximum value for y axis
#' @param tick_dist Distance in basepairs between chromosome diagram ticks
#' @param axis_label_margin Size of margin between axis and its label (distance between axis and axis label)
#' @param main_title_cex Scaling factor for main title
#' @param axis_label_cex Scaling factor for axis label
#' @param axis_label_offset Offset for axis label
#' @param axis_tick_cex Scaling factor for axis ticks
#' @param plot_genes Flag to have karyoploteR plot gene structures beneath Manhattan
#' @param gene_names_cex Scaling factor for gene labels
#' @param chrom_names_cex Scaling factor for chromosome labels
#' @param third_gRanges GenomicRanges object containing data for third plot to be placed above back-to-back Manhattans
#' @param third_label Axis label for third plot
#' @param chromosomes List of chromosomes to plot. Expects strings of the form 'chrx' where 'x' is from the set {1, 2, ..., 22, X, Y} 
#' @param zoom Coordinates for interval to magnify. Expects a string of the form 'chrx:a-b' where x is the chromosome number a and b are basepair coordinates.
#' @param plot_params List of named plot parameters to pass to \code{plotKaryotype}
#' @param top_gRanges_points.col Colours used to plot the points in top_gRanges
#' @param bottom_gRanges_points.col Colours used to plot the points in bottom_gRanges
#' @param third_gRanges_points.col Colours used to plot the points in third_gRanges
#' @param points.cex Size of the point symbols
#' @param top_highlight (GRanges, character vector, logical vector or numeric vector) The points to highlight in a different color in the top plot. If a GRanges (or anythng accepted by toGRanges) the points overlapping these regions will be highlighted. Otherwise the points will be selected with data[highlight]. If NULL no point will be highlighted. (defaults to NULL)
#' @param bottom_highlight (GRanges, character vector, logical vector or numeric vector) The points to highlight in a different color in the bottom plot. If a GRanges (or anythng accepted by toGRanges) the points overlapping these regions will be highlighted. Otherwise the points will be selected with data[highlight]. If NULL no point will be highlighted. (defaults to NULL)
#' @param third_highlight (GRanges, character vector, logical vector or numeric vector) The points to highlight in a different color in the third plot. If a GRanges (or anythng accepted by toGRanges) the points overlapping these regions will be highlighted. Otherwise the points will be selected with data[highlight]. If NULL no point will be highlighted. (defaults to NULL)
#' @param top_highlight.col Top plot highlight colour
#' @param bottom_highlight.col Bottom plot highlight colour
#' @param third_highlight.col Third plot highlight colour
#'
#' @return KaryoPlot object
#'
#' @importFrom karyoploteR plotKaryotype kpAddBaseNumbers kpAddChromosomeNames kpAddLabels kpAxis kpPlotGenes kpPlotManhattan makeGenesDataFromTxDb addGeneNames mergeTranscripts getDefaultPlotParams
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom GenomicRanges mcols mcols<-
#' @export
back_to_back_manhattan <- function(top_gRanges, bottom_gRanges, top_label, bottom_label, main, ymax = 15, tick_dist = 1e5, axis_label_margin = 0.03, main_title_cex = 2.7, axis_label_cex = 1.8, axis_label_offset = 0, axis_tick_cex = 1, plot_genes = F, gene_names_cex = 1, chrom_names_cex = 2, third_gRanges = NULL, third_label = NULL, chromosomes = NULL, zoom = NULL, plot_params = getDefaultPlotParams(plot.type = 4), top_gRanges_points.col = '2blues', bottom_gRanges_points.col = '2blues', third_gRanges_points.col = '2blues', points.cex = 1, top_highlight = NULL, bottom_highlight = NULL, third_highlight = NULL, top_highlight.col = 'greenyellow', bottom_highlight.col = 'greenyellow', third_highlight.col = 'greenyellow') {

  if(!is.null(chromosomes) & !is.null(zoom)) {
    stop("Cannot specify both \'chromosomes\' and \'zoom\'")
  }

  if(!is.null(third_gRanges)) {
    if(is.null(third_label)) {
      stop("Need to specify \'third_label\' if specifying \'third_gRanges\'")
      }
    if(!('P' %in% names(mcols(third_gRanges)) | 'p' %in% names(mcols(third_gRanges)))) {
      stop('Need a p-value column labelled \'P\' or \'p\' in third_gRanges')
    }
  }

  if(plot_genes & is.null(zoom)) {
    stop("Can only plot gene tracks if \'zoom\' coordinates are specified")
  }

  if(!('P' %in% names(mcols(top_gRanges)) | 'p' %in% names(mcols(top_gRanges)))) {
    stop('Need a p-value column labelled \'P\' or \'p\' in top_gRanges')
  }

  if(!('P' %in% names(mcols(bottom_gRanges)) | 'p' %in% names(mcols(bottom_gRanges)))) {
    stop('Need a p-value column labelled \'P\' or \'p\' in bottom_gRanges')
  } 

  if(!is.null(chromosomes)) {
    kp <- plotKaryotype(plot.type = 4, labels.plotter = NULL, chromosomes = chromosomes, plot.params = plot_params)
    kp <- kpAddBaseNumbers(kp, add.units = T, cex = 1, tick.dist = tick_dist)
  } else if(!is.null(zoom)) {
    kp <- plotKaryotype(plot.type = 4, labels.plotter = NULL, zoom = zoom, plot.params = plot_params)
    kp <- kpAddBaseNumbers(kp, add.units = T, cex = 1, tick.dist = tick_dist)
  } else {
    kp <- plotKaryotype(plot.type = 4, labels.plotter = NULL, plot.params = plot_params)
    kp <- kpAddChromosomeNames(kp, col = 'black',srt = 90,cex = chrom_names_cex)
  } 

  title(main = main, cex.main =  main_title_cex)

  if(is.null(third_gRanges)) {
    if(plot_genes) {
      # two tracks, plot genes
      kp <- kpAddLabels(kp, labels = top_label, srt=90, pos=3, r0=0.6, r1=1, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.6, r1=1, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data = top_gRanges, r0 = 0.6, r1 = 1, ymax = ymax, points.cex = points.cex, points.col = top_gRanges_points.col, highlight = top_highlight, highlight.col = top_highlight.col)

      kp <- kpAddLabels(kp, labels = bottom_label, srt=90, pos=3, r0=0.6, r1=0.2, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.6, r1=0.2, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=bottom_gRanges, r0=0.6, r1=0.2, ymax=ymax, points.cex = points.cex, points.col = bottom_gRanges_points.col, highlight = bottom_highlight, highlight.col = bottom_highlight.col)

    } else {
      # two tracks, do not plot genes
      kp <- kpAddLabels(kp, labels = top_label, srt=90, pos=3, r0=0.5, r1=1, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.5, r1=1, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=top_gRanges, r0=0.5, r1=1, ymax=ymax, points.cex = points.cex, points.col = top_gRanges_points.col, highlight = top_highlight, highlight.col = top_highlight.col)


      kp <- kpAddLabels(kp, labels = bottom_label, srt=90, pos=3, r0=0, r1=0.5, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.5, r1=0, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=bottom_gRanges, r0=0.5, r1=0, ymax=ymax, points.cex = points.cex, points.col = bottom_gRanges_points.col, highlight = bottom_highlight, highlight.col = bottom_highlight.col)
    }
  } else {
    if(plot_genes) {
      # Three-track Manhattan, plot genes
      kp <- kpAddLabels(kp, labels = third_label, srt=90, pos=3, r0=0.8, r1=1, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.8, r1=1, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=third_gRanges, r0=0.8, r1=1, ymax=ymax, points.cex = points.cex, points.col = third_gRanges_points.col, highlight = third_highlight, highlight.col = third_highlight.col)


      kp <- kpAddLabels(kp, labels = top_label, srt=90, pos=3, r0=0.46, r1=0.72, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.46, r1=0.72, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=top_gRanges, r0=0.46, r1=0.72, ymax=ymax, points.cex = points.cex, points.col = top_gRanges_points.col, highlight = top_highlight, highlight.col = top_highlight.col)


      kp <- kpAddLabels(kp, labels = bottom_label, srt=90, pos=3, r0=0.2, r1=0.46, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.46, r1=0.2, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=bottom_gRanges, r0=0.46, r1=0.2, ymax=ymax, points.cex = points.cex, points.col = bottom_gRanges_points.col, highlight = bottom_highlight, highlight.col = bottom_highlight.col)

    } else {
      # Three-track Manhattan
      kp <- kpAddLabels(kp, labels = third_label, srt=90, pos=3, r0=0.7, r1=1, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.7, r1=1, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=third_gRanges, r0=0.7, r1=1, ymax=ymax, points.cex = points.cex, points.col = third_gRanges_points.col, highlight = third_highlight, highlight.col = third_highlight.col)

      kp <- kpAddLabels(kp, labels = top_label, srt=90, pos=3, r0=0.33, r1=0.66, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.33, r1=0.66, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=top_gRanges, r0=0.33, r1=0.66, ymax=ymax, points.cex = points.cex, points.col = top_gRanges_points.col, highlight = top_highlight, highlight.col = top_highlight.col)

      kp <- kpAddLabels(kp, labels = bottom_label, srt=90, pos=3, r0=0, r1=0.33, cex=axis_label_cex, label.margin = axis_label_margin, offset = axis_label_offset)
      # Note how r0 and r1 are flipped here for kpAxis and kpPlotManhattan
      kp <- kpAxis(kp, ymin=0, ymax=ymax, r0=0.33, r1=0, cex = axis_tick_cex)
      kp <- kpPlotManhattan(kp, data=bottom_gRanges, r0=0.33, r1=0, ymax=ymax, points.cex = points.cex, points.col = bottom_gRanges_points.col, highlight = bottom_highlight, highlight.col = bottom_highlight.col)
      }
  }

  if(plot_genes) {
    genes.data <- makeGenesDataFromTxDb(karyoplot=kp, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)

    genes.data <- addGeneNames(genes.data)
    genes.data.merged <- mergeTranscripts(genes.data)
    kp <- kpPlotGenes(kp, data=genes.data.merged, r0=0, r1=0.18, gene.name.cex=gene_names_cex, gene.name.position='left')
  }

  kp
}
