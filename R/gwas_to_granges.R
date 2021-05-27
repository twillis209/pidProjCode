#' @title Convert GWAS data in data.table to GenomicRanges object
#' 
#' @param data_table data.table containing GWAS data
#' @param bp_label name of basepair column
#' @param chr_label name of chromosome column
#' @param max_p upper bound on p-values
#' @param min_p lower bound on p-values
#' 
#' @return GenomicRanges object
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
#'
#' TODO upper bound is not inclusive, i.e. <= does not seem to work
#' @export
gwas_to_granges <- function(data_table, bp_label = 'BP38', chr_label = 'CHR38', max_p = NULL, min_p = NULL) {
  chr_temp <- data_table[[chr_label]]

  data_table[, (chr_label) := paste0('chr', get(chr_label))]

  granges <- makeGRangesFromDataFrame(data.frame(data_table), start.field = bp_label, end.field = bp_label, seqnames.field = chr_label, ignore.strand = T, keep.extra.columns = T)

  data_table[, (chr_label) := chr_temp]

  if(!is.null(max_p)) {
    if('P' %in% names(mcols(granges))) {
      granges <- granges[mcols(granges)$P <= max_p]
    } else if('p' %in% names(mcols(granges))) {
      granges <- granges[mcols(granges)$p <= max_p]
    } else {
      stop("max_p specified but neither \'P\' nor \'p\' are names of metadata columns")
    }
  }

  if(!is.null(min_p)) {
    if('P' %in% names(mcols(granges))) {
      granges <- granges[mcols(granges)$P >= min_p]
    } else if('p' %in% names(mcols(granges))) {
      granges <- granges[mcols(granges)$p >= min_p]
    } else {
      stop("min_p specified but neither \'P\' nor \'p\' are names of metadata columns")
    }
  }

  granges
}
