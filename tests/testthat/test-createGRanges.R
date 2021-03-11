test_that("GRanges object is created in simple use case", {
  set.seed(42)
  gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  granges<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp')
})

test_that("GRanges object is created with correct names", {
  set.seed(42)
  gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  granges<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp', nameCol='SNP')
  # GenomicRanges sorts the names in an odd way
  expect_equal(names(granges), sort(paste0('rs',1:10))) 
})

test_that("Error is thrown for chrCol argument not found in dataFrame", {
  set.seed(42)
  gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  expect_error(createGRanges(gwasDaf, chrCol='chr1', bpCol='bp'), 'chrCol name is not a column name in dataFrame')
})

test_that("Error is thrown for bpCol argument not found in dataFrame", {
  set.seed(42)
  gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  expect_error(createGRanges(gwasDaf, chrCol='chr', bpCol='bp1'), 'bpCol name is not a column name in dataFrame')
})

test_that("Error is thrown for nameCol argument not found in dataFrame", {
  set.seed(42)
  gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  expect_error(createGRanges(gwasDaf, chrCol='chr', bpCol='bp', nameCol='missing'), 'nameCol name is not a column name in dataFrame')
})
