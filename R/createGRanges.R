# Code is adapted from cpvSNP::createArrayData
#' @title Create a GenomicRanges::GRanges object.
#'
#' @description Creates a GenomicRanges::GRanges object from a data frame containing genomic coordinates and metadata.
#'
#' @details Adapted from the code for cpvSNP::createArrayData.
#'
#' @param dataFrame A data.frame containing genomic coordinates and metadata
#' @param chrCol A string giving the name of chromosome column
#' @param bpCol A string giving the name of basepair column
#'
#' @return A GenomicRanges::GRanges object containing the same (meta)data as dataFrame
#' 
#' @importFrom GenomicRanges GRanges elementMetadata<- 
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
#' @importFrom methods is
#' @export
#'
#' @examples
#' 
#' gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
#' granges<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp')
#' 
createGRanges<-function(dataFrame, chrCol, bpCol) {
  if(!is(dataFrame, 'data.frame') | is(dataFrame, 'data.table')) {
    stop('dataFrame must be a data.frame object')
  }

  if(!any(is.element(names(dataFrame), chrCol))) {
    stop('chrCol name is not a column name in dataFrame')
  }

  if(!any(is.element(names(dataFrame), bpCol))) {
    stop('bpCol name is not a column name in dataFrame')
  }

  names(dataFrame)[is.element(names(dataFrame), chrCol)]<-'chromosome'

  dataFrame$Start<-dataFrame[,bpCol]
  dataFrame$End<-dataFrame[,bpCol]

  dataFrame$Start[is.na(dataFrame$chromosome)]<-1
  dataFrame$End[is.na(dataFrame$chromosome)]<-1

  dataFrame$chromosome[is.na(dataFrame$chromosome)]<-'U'

  dataFrame$chromosome <- factor(dataFrame$chromosome)

  levels(dataFrame$chromosome)<-paste('chr', levels(dataFrame$chromosome), sep = '')

  dataFrame<-dataFrame[order(dataFrame$chromosome, dataFrame$Start),]

  granges<-GRanges(Rle(dataFrame$chromosome),IRanges(start=dataFrame$Start, end=dataFrame$End))

  dropIndices<-which(names(dataFrame) %in% c('chromosome','Start','End',bpCol))

  elementMetadata(granges)<-dataFrame[,-dropIndices]

  granges
}
