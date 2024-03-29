#' @title Writes out a multitrack Manhattan.
#' 
#' @description Writes out a list of \code{GenomicRanges} objects in a stacked Manhattan plot.
#'
#'
#' @param gRanges List of \code{GenomicRanges} objects to plot
#' @param axis_labels List of axis labels
#' @param main Main title
#' @param ymax Maximum value for y axes
#' @param axis_label_margin Size of margin between axis and its label (distance between axis and axis label)
#' @param main_title_cex Scaling factor for main title
#' @param axis_label_cex Scaling factor for axis label
#' @param axis_label_offset Offset for axis label
#' @param axis_tick_cex Scaling factor for axis ticks
#' @param chrom_names_cex Scaling factor for chromosome labels
#' @param chromosomes List of chromosomes to plot. Expects strings of the form 'chrx' where 'x' is from the set {1, 2, ..., 22, X, Y} 
#' @param chrom_tick_dist Distance in bp between ideogram ticks
#' @param chrom_tick_cex Scaling factor for ideogram ticks
#' @param plot_params List of named plot parameters to pass to \code{plotKaryotype}
#' @param track_margin Proportion of track assigned to track margin
#' @param points.col Colours used to plot points
#' @param points.cex Size of the point symbols
#' @param highlights List of GRanges, or character, numeric, or logical vectors giving the points to highlight in a different color in the top plot. If a GRanges (or anythng accepted by toGRanges) the points overlapping these regions will be highlighted. Otherwise the points will be selected with data[highlight]. If NULL no point will be highlighted. (defaults to NULL) 
#' @param highlight.cols Colours for each highlight
#'
#' @return KaryoPlot
#'
#' @importFrom karyoploteR plotKaryotype kpAddBaseNumbers kpAddChromosomeNames kpAddLabels kpAxis kpPlotManhattan getDefaultPlotParams autotrack
#' @importFrom GenomicRanges mcols mcols<-
#' @export
multitrack_manhattan <- function(gRanges, axis_labels, main, ymax = 15, axis_label_margin = 0.03, main_title_cex = 2.7, axis_label_cex = 1.8, axis_label_offset = 0, axis_tick_cex = 1, chrom_names_cex = 2, chromosomes = NULL, chrom_tick_dist=1e6, chrom_tick_cex = 1, plot_params = getDefaultPlotParams(plot.type = 4), track_margin = 0.06, points.col = '2blues', points.cex = 1, highlights = NULL, highlight.cols = NULL) {

  for(i in seq_along(gRanges)) {
    if(!any(c('P', 'p') %in% names(mcols(gRanges[[i]])))) {
      stop(sprintf("Element %d of gRanges does not contain a column labelled \'P\' or \'p\'", i))
      }
    }

  if(!is.null(chromosomes)) {
    kp <- plotKaryotype(plot.type = 4, labels.plotter = NULL, chromosomes = chromosomes, plot.params = plot_params)
    kp <- kpAddBaseNumbers(kp, add.units = T, cex = chrom_tick_cex, tick.dist = chrom_tick_dist)
  } else {
    kp <- plotKaryotype(plot.type = 4, labels.plotter = NULL, plot.params = plot_params, chromosomes = paste0('chr', 1:22))
    kp <- kpAddChromosomeNames(kp, col = 'black', srt = 90, cex = chrom_names_cex)
  }

  title(main = main, cex.main = main_title_cex)

  for(i in seq_along(gRanges)) {
    auto <- autotrack(current.track = i, total.tracks = length(gRanges), margin = track_margin)

    kp <- kpAddLabels(kp, labels = axis_labels[i],
                srt = 90, pos = 3,
                r0 = auto$r0, r1 = auto$r1,
                cex = axis_label_cex,
                label.margin = axis_label_margin,
                offset = axis_label_offset)

    kp <- kpAxis(kp, ymin = 0, ymax = ymax,
           r0 = auto$r0, r1 = auto$r1, cex = axis_tick_cex)
    if(is.null(highlights)) {
      kp <- kpPlotManhattan(kp, data = gRanges[[i]],
                            points.cex = points.cex,
                          points.col = points.col,
                          r0 = auto$r0, r1 = auto$r1,
                          ymax = ymax)
    } else {
      kp <- kpPlotManhattan(kp, data = gRanges[[i]],
                            points.cex = points.cex,
                            points.col = points.col,
                            r0 = auto$r0, r1 = auto$r1,
                            ymax = ymax,
                            highlight = highlights[[i]],
                            highlight.col = highlight.cols[[i]])
      }
  }

  kp
}
