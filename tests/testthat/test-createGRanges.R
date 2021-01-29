test_that("GRanges object is created in simple use case", {
  set.seed(42)
  gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  granges<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp')
})
