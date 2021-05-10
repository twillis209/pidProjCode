#' @title Add names of genes local to ranges stored in GenomicRanges object
#'
#' @param granges GenomicRanges object
#' @param txdb TxDb containing annotations
#' @param flank_length width of flanking region to add to each range
#'
#' @return GenomicRanges with gene names added
#'
#' @export
add_gene_names_to_granges <- function(granges, txdb, flank_length = NULL) {
  granges$home_gene <- get_gene_names_for_granges(granges, txdb)

  if(!is.null(flank_length)) {
    granges$genes <- get_gene_names_for_granges(granges, txdb, flank_length = flank_length)
  }

  granges
}

#' @title Fetch genes from TxDb for GenomicRanges object
#'
#' @details Adapted from http://www.bioinsteps.com/2019/08/annotate-genomic-coordinates-using-r.html
#'
#' @param granges GenomicRanges object
#' @param txdb TxDb containing annotations
#' @param flank_length width of flanking region to add to each range
#'
#' @importFrom GenomicFeatures genes
#' @importFrom annotate getSYMBOL
#' @importFrom IRanges findOverlaps flank
#' @importFrom S4Vectors queryHits queryLength splitAsList subjectHits
#' @importFrom org.Hs.eg.db org.Hs.eg
#'
#' @export
get_gene_names_for_granges <- function(granges, txdb, flank_length = NULL) {

  genes <- genes(txdb)

  if(!is.null(flank_length)) {
    granges <- flank(granges, flank_length, both = T)
  }

  olaps <- findOverlaps(granges, genes)
  mcols(olaps)$gene_id <- genes$gene_id[subjectHits(olaps)]
  granges_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
  gene_id <- splitAsList(mcols(olaps)$gene_id, granges_factor)

  gene_names <- list()

  for(i in 1:length(granges)) {
    if(length(na.omit(gene_id[[i]])) != 0) {
      # TODO does this need an import?
      gene_names[[i]] <- unname(getSYMBOL(gene_id[[i]], data = 'org.Hs.eg'))
    }
  }

  gene_names
}

#' @title Fetch CDS ID from TxDb for GenomicRanges object
#'
#' @param granges GenomicRanges object
#' @param txdb TxDb containing annotations
#'
#' @importFrom GenomicFeatures cds
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits queryLength splitAsList subjectHits
#'
#' @export
get_cds_id_for_granges <- function(granges, txdb) {
  cds_granges <- cds(txdb)

  olaps <- findOverlaps(granges, cds_granges)
  mcols(olaps)$cds_id <- cds_granges$cds_id[subjectHits(olaps)]
  granges_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))

  splitAsList(mcols(olaps)$cds_id, granges_factor)
}

#' @title Fetch exon ID from TxDb for GenomicRanges object
#'
#' @param granges GenomicRanges object
#' @param txdb TxDb containing annotations
#'
#' @importFrom GenomicFeatures exons
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits queryLength splitAsList subjectHits
#'
#' @export
get_exon_id_for_granges <- function(granges, txdb) {
  exon_granges <- exons(txdb)

  olaps <- findOverlaps(granges, exon_granges)
  mcols(olaps)$exon_id <- exon_granges$exon_id[subjectHits(olaps)]
  granges_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))

  splitAsList(mcols(olaps)$exon_id, granges_factor)
}
