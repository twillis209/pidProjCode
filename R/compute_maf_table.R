#' @title Returns table of MAF values by interval.
#'
#' @param maf Vector of MAF values
#' @param breaks Vector of breaks to define intervals for table
#' @param proportions Flag to return table with proportions rather than counts
#'
#' @return Table
#' @export
#'
compute_maf_table <- function(maf, breaks = seq(0, 0.5, length = 51), proportions = F) {
  maf_interval <- addNA(cut(maf, breaks = breaks, include.lowest = T))

  maf_interval_freq_whole <- table(maf_interval)

  if(proportions) {
    maf_interval_freq_whole / sum(maf_interval_freq_whole)
  } else {
    maf_interval_freq_whole
  }
}
