test_that("Error is thrown for chrCol argument not found in dataFrame", {
  set.seed(42)
  gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  grangesTop<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp')
  grangesBottom<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp')
  expect_error(backToBackManhattan(topGRanges=grangesTop,
                                   bottomGRanges=grangesBottom,
                                   outputFile='test.png',
                                   topLabel='top',
                                   bottomLabel='bottom',
                                   main='title',
                                   chromosomes=c('chr1'),
                                   zoom='chr1:28567-28576'),
               'Cannot specify both \'chromosomes\' and \'zoom\'')
})

test_that("Error is thrown for thirdGRanges argument without accompanying thirdLabel argument", {
  set.seed(42)
  gwasDaf<-data.frame('chr'=1:10, 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  grangesTop<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp')
  grangesBottom<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp')
  grangesThird<-createGRanges(gwasDaf, chrCol='chr', bpCol='bp')
  expect_error(backToBackManhattan(topGRanges=grangesTop,
                                   bottomGRanges=grangesBottom,
                                   outputFile='test.png',
                                   topLabel='top',
                                   bottomLabel='bottom',
                                   main='title',
                                   thirdGRanges=grangesThird,
                                   thirdLabel=NULL),
               'Need to specify \'thirdLabel\' if specifying \'thirdGRanges\'')
})
