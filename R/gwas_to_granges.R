#' @title Convert GWAS data in data.table to GenomicRanges object
#' 
#' @param data_table data.table containing GWAS data
#' @param bp_label name of basepair column
#' @param chr_label name of chromosome column
#' 
#' @return GenomicRanges object
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
#' 
#' @export
gwas_to_granges <- function(data_table, bp_label = 'BP38', chr_label = 'CHR38') {
  data_table[, (chr_label) := paste0('chr', get(chr_label))]

  makeGRangesFromDataFrame(data.frame(data_table), start.field = bp_label, end.field = bp_label, seqnames.field = chr_label, ignore.strand = T, keep.extra.columns = T)
}
