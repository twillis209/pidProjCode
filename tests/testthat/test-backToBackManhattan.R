test_that("Error is thrown when both \'chromosomes\' and \'zoom\' arguments are specified", {
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

test_that("Error is thrown for specifying \'plotGenes\' flag without accompanying \'zoom\' argument", {
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
                                   plotGenes=T),
               'Can only plot gene tracks if \'zoom\' coordinates are specified')
})

test_that("Whole-genome Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidGRanges <- createGRanges(data.frame(pidDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  igadGRanges <- createGRanges(data.frame(igadDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')

  backToBackManhattan(pidGRanges, igadGRanges, output='test_all.png', topLabel='PID', bottomLabel='IgAD', main='PID and IgAD')

  expect_equal(2, 1+1)
})

test_that("Whole-karyotype Manhattan is drawn when outputFile is NULL", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidGRanges <- createGRanges(data.frame(pidDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  igadGRanges <- createGRanges(data.frame(igadDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')

  backToBackManhattan(pidGRanges, igadGRanges, outputFile=NULL, topLabel='PID', bottomLabel='IgAD', main='PID and IgAD')

  expect_equal(2, 1+1)
})

test_that("Single-chromosome Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidGRanges <- createGRanges(data.frame(pidDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  igadGRanges <- createGRanges(data.frame(igadDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')

  backToBackManhattan(pidGRanges, igadGRanges, output='test_chr6.png', topLabel='PID', bottomLabel='IgAD', main='Chromosome 6', chromosomes='chr6', tickDist=1e7)

  expect_equal(2, 1+1)
})

test_that("Single-chromosome Manhattan with three tracks is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  sleDat <- data.table::fread(file=system.file('extdata', 'sleGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(sleDat, -log10(P) < 15)

  pidGRanges <- createGRanges(data.frame(pidDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  igadGRanges <- createGRanges(data.frame(igadDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  sleGRanges <- createGRanges(data.frame(sleDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')

  backToBackManhattan(pidGRanges, igadGRanges, output='test_chr6_three_track.png', topLabel='PID', bottomLabel='IgAD', thirdLabel = 'SLE', thirdGRanges = sleGRanges, main='Chromosome 6', chromosomes='chr6', tickDist=1e7)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidGRanges <- createGRanges(data.frame(pidDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  igadGRanges <- createGRanges(data.frame(igadDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')

  backToBackManhattan(pidGRanges, igadGRanges, output='test_chr6_zoom.png', topLabel='PID', bottomLabel='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-36e6', tickDist=1e6)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan is drawn with genes", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidGRanges <- createGRanges(data.frame(pidDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  igadGRanges <- createGRanges(data.frame(igadDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')

  backToBackManhattan(pidGRanges, igadGRanges, output='test_chr6_zoom_genes.png', topLabel='PID', bottomLabel='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-26e6', tickDist=1e5, plotGenes=T)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan with three plots is drawn with genes", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  sleDat <- data.table::fread(file=system.file('extdata', 'sleGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)
  sleDat <- subset(sleDat, -log10(P) < 15)

  pidGRanges <- createGRanges(data.frame(pidDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  igadGRanges <- createGRanges(data.frame(igadDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')
  sleGRanges <- createGRanges(data.frame(sleDat), chrCol='CHR38', bpCol='BP38', nameCol='SNPID')

  backToBackManhattan(pidGRanges, igadGRanges, thirdGRanges=sleGRanges, thirdLabel='SLE', output='test_chr6_zoom_genes_three_track.png', topLabel='PID', bottomLabel='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-26e6', tickDist=1e5, plotGenes=T)

  expect_equal(2, 1+1)
})
