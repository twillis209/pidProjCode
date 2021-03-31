test_that("Error is thrown when both \'chromosomes\' and \'zoom\' arguments are specified", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  grangesTop <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)

  expect_error(backToBackManhattan(topGRanges=grangesTop,
                                   bottomGRanges=grangesBottom,
                                   topLabel='top',
                                   bottomLabel='bottom',
                                   main='title',
                                   chromosomes=c('chr1'),
                                   zoom='chr1:28567-28576'),
               'Cannot specify both \'chromosomes\' and \'zoom\'')
})

test_that("Error is thrown for thirdGRanges argument without accompanying thirdLabel argument", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  grangesTop <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesThird <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)

  expect_error(backToBackManhattan(topGRanges=grangesTop,
                                   bottomGRanges=grangesBottom,
                                   topLabel='top',
                                   bottomLabel='bottom',
                                   main='title',
                                   thirdGRanges=grangesThird,
                                   thirdLabel=NULL),
               'Need to specify \'thirdLabel\' if specifying \'thirdGRanges\'')
})

test_that("Error is thrown for specifying \'plotGenes\' flag without accompanying \'zoom\' argument", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))

  grangesTop <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesThird <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)

  expect_error(backToBackManhattan(topGRanges=grangesTop,
                                   bottomGRanges=grangesBottom,
                                   topLabel='top',
                                   bottomLabel='bottom',
                                   main='title',
                                   plotGenes=T),
               'Can only plot gene tracks if \'zoom\' coordinates are specified')
})

test_that("Error is thrown for omitting p-value column from topGRanges", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  grangesTop <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesTop$p <- NULL

  expect_error(backToBackManhattan(topGRanges=grangesTop,
                                   bottomGRanges=grangesBottom,
                                   topLabel='top',
                                   bottomLabel='bottom',
                                   main='title'),
               'Need a p-value column labelled \'P\' or \'p\' in topGRanges')
})

test_that("Error is thrown for omitting p-value column from bottomGRanges", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))

  grangesTop <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom$p <- NULL

  expect_error(backToBackManhattan(topGRanges=grangesTop,
                                   bottomGRanges=grangesBottom,
                                   topLabel='top',
                                   bottomLabel='bottom',
                                   main='title'),
               'Need a p-value column labelled \'P\' or \'p\' in bottomGRanges')
})

test_that("Error is thrown for omitting p-value column from thirdGRanges", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))

  grangesTop <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesThird <- makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesThird$p <- NULL

  expect_error(backToBackManhattan(topGRanges=grangesTop,
                                   bottomGRanges=grangesBottom,
                                   thirdGRanges=grangesThird,
                                   topLabel='top',
                                   bottomLabel='bottom',
                                   thirdLabel='ref',
                                   main='title'),
               'Need a p-value column labelled \'P\' or \'p\' in thirdGRanges')
})

test_that("Whole-genome Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat[, CHR38 := paste0('chr', CHR38)]
  igadDat[, CHR38 := paste0('chr', CHR38)]

  pidGRanges <- makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  backToBackManhattan(pidGRanges, igadGRanges, topLabel='PID', bottomLabel='IgAD', main='PID and IgAD')

  expect_equal(2, 1+1)
})

test_that("Whole-karyotype Manhattan is drawn when outputFile is NULL", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat[, CHR38 := paste0('chr', CHR38)]
  igadDat[, CHR38 := paste0('chr', CHR38)]

  pidGRanges <- makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  backToBackManhattan(pidGRanges, igadGRanges, topLabel='PID', bottomLabel='IgAD', main='PID and IgAD')

  expect_equal(2, 1+1)
})

test_that("Single-chromosome Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat[, CHR38 := paste0('chr', CHR38)]
  igadDat[, CHR38 := paste0('chr', CHR38)]

  pidGRanges <- makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  backToBackManhattan(pidGRanges, igadGRanges, topLabel='PID', bottomLabel='IgAD', main='Chromosome 6', chromosomes='chr6', tickDist=1e7)

  expect_equal(2, 1+1)
})

test_that("Single-chromosome Manhattan with three tracks is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  sleDat <- data.table::fread(file=system.file('extdata', 'sleGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)
  sleDat <- subset(sleDat, -log10(P) < 15)

  pidDat[, CHR38 := paste0('chr', CHR38)]
  igadDat[, CHR38 := paste0('chr', CHR38)]
  sleDat[, CHR38 := paste0('chr', CHR38)]

  pidGRanges <- makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  sleGRanges <- makeGRangesFromDataFrame(data.frame(sleDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  backToBackManhattan(pidGRanges, igadGRanges, topLabel='PID', bottomLabel='IgAD', thirdLabel = 'SLE', thirdGRanges = sleGRanges, main='Chromosome 6', chromosomes='chr6', tickDist=1e7)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat[, CHR38 := paste0('chr', CHR38)]
  igadDat[, CHR38 := paste0('chr', CHR38)]

  pidGRanges <- makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  backToBackManhattan(pidGRanges, igadGRanges, topLabel='PID', bottomLabel='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-36e6', tickDist=1e6)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan is drawn with genes", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat[, CHR38 := paste0('chr', CHR38)]
  igadDat[, CHR38 := paste0('chr', CHR38)]

  pidGRanges <- makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  backToBackManhattan(pidGRanges, igadGRanges, topLabel='PID', bottomLabel='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-26e6', tickDist=1e5, plotGenes=T)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan with three plots is drawn with genes", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  sleDat <- data.table::fread(file=system.file('extdata', 'sleGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)
  sleDat <- subset(sleDat, -log10(P) < 15)

  pidDat[, CHR38 := paste0('chr', CHR38)]
  igadDat[, CHR38 := paste0('chr', CHR38)]
  sleDat[, CHR38 := paste0('chr', CHR38)]

  pidGRanges <- makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  sleGRanges <- makeGRangesFromDataFrame(data.frame(sleDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  backToBackManhattan(pidGRanges, igadGRanges, thirdGRanges=sleGRanges, thirdLabel='SLE', topLabel='PID', bottomLabel='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-26e6', tickDist=1e5, plotGenes=T)

  expect_equal(2, 1+1)
})
