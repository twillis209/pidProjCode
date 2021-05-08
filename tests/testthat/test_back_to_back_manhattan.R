test_that("Error is thrown when both \'chromosomes\' and \'zoom\' arguments are specified", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  grangesTop <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)

  expect_error(back_to_back_manhattan(top_gRanges=grangesTop,
                                   bottom_gRanges=grangesBottom,
                                   top_label='top',
                                   bottom_label='bottom',
                                   main='title',
                                   chromosomes=c('chr1'),
                                   zoom='chr1:28567-28576'),
               'Cannot specify both \'chromosomes\' and \'zoom\'')
})

test_that("Error is thrown for third_gRanges argument without accompanying third_label argument", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  grangesTop <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesThird <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)

  expect_error(back_to_back_manhattan(top_gRanges=grangesTop,
                                   bottom_gRanges=grangesBottom,
                                   top_label='top',
                                   bottom_label='bottom',
                                   main='title',
                                   third_gRanges=grangesThird,
                                   third_label=NULL),
               'Need to specify \'third_label\' if specifying \'third_gRanges\'')
})

test_that("Error is thrown for specifying \'plot_genes\' flag without accompanying \'zoom\' argument", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))

  grangesTop <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesThird <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)

  expect_error(back_to_back_manhattan(top_gRanges=grangesTop,
                                   bottom_gRanges=grangesBottom,
                                   top_label='top',
                                   bottom_label='bottom',
                                   main='title',
                                   plot_genes=T),
               'Can only plot gene tracks if \'zoom\' coordinates are specified')
})

test_that("Error is thrown for omitting p-value column from top_gRanges", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))
  grangesTop <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesTop$p <- NULL

  expect_error(back_to_back_manhattan(top_gRanges=grangesTop,
                                   bottom_gRanges=grangesBottom,
                                   top_label='top',
                                   bottom_label='bottom',
                                   main='title'),
               'Need a p-value column labelled \'P\' or \'p\' in top_gRanges')
})

test_that("Error is thrown for omitting p-value column from bottom_gRanges", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))

  grangesTop <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom$p <- NULL

  expect_error(back_to_back_manhattan(top_gRanges=grangesTop,
                                   bottom_gRanges=grangesBottom,
                                   top_label='top',
                                   bottom_label='bottom',
                                   main='title'),
               'Need a p-value column labelled \'P\' or \'p\' in bottom_gRanges')
})

test_that("Error is thrown for omitting p-value column from third_gRanges", {
  set.seed(42)
  gwasDaf <- data.frame('chr'=paste0('chr', 1:10), 'bp'=28567:28576, 'SNP'=paste0('rs', 1:10),'p'=runif(10))

  grangesTop <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesBottom <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesThird <- GenomicRanges::makeGRangesFromDataFrame(gwasDaf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  grangesThird$p <- NULL

  expect_error(back_to_back_manhattan(top_gRanges=grangesTop,
                                   bottom_gRanges=grangesBottom,
                                   third_gRanges=grangesThird,
                                   top_label='top',
                                   bottom_label='bottom',
                                   third_label='ref',
                                   main='title'),
               'Need a p-value column labelled \'P\' or \'p\' in third_gRanges')
})

test_that("Whole-genome Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat$CHR38 <- paste0('chr', pidDat$CHR38)
  igadDat$CHR38 <- paste0('chr', igadDat$CHR38)

  pidGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  back_to_back_manhattan(pidGRanges, igadGRanges, top_label='PID', bottom_label='IgAD', main='PID and IgAD')

  expect_equal(2, 1+1)
})

test_that("Whole-karyotype Manhattan is drawn when outputFile is NULL", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat$CHR38 <- paste0('chr', pidDat$CHR38)
  igadDat$CHR38 <- paste0('chr', igadDat$CHR38)

  pidGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  back_to_back_manhattan(pidGRanges, igadGRanges, top_label='PID', bottom_label='IgAD', main='PID and IgAD')

  expect_equal(2, 1+1)
})

test_that("Single-chromosome Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat$CHR38 <- paste0('chr', pidDat$CHR38)
  igadDat$CHR38 <- paste0('chr', igadDat$CHR38)

  pidGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  back_to_back_manhattan(pidGRanges, igadGRanges, top_label='PID', bottom_label='IgAD', main='Chromosome 6', chromosomes='chr6', tick_dist=1e7)

  expect_equal(2, 1+1)
})

test_that("Single-chromosome Manhattan with three tracks is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  sleDat <- data.table::fread(file=system.file('extdata', 'sleGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)
  sleDat <- subset(sleDat, -log10(P) < 15)

  pidDat$CHR38 <- paste0('chr', pidDat$CHR38)
  igadDat$CHR38 <- paste0('chr', igadDat$CHR38)
  sleDat$CHR38 <- paste0('chr', sleDat$CHR38)

  pidGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  sleGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(sleDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  back_to_back_manhattan(pidGRanges, igadGRanges, top_label='PID', bottom_label='IgAD', third_label = 'SLE', third_gRanges = sleGRanges, main='Chromosome 6', chromosomes='chr6', tick_dist=1e7)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan is drawn", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat$CHR38 <- paste0('chr', pidDat$CHR38)
  igadDat$CHR38 <- paste0('chr', igadDat$CHR38)

  pidGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  back_to_back_manhattan(pidGRanges, igadGRanges, top_label='PID', bottom_label='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-36e6', tick_dist=1e6)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan is drawn with genes", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)

  pidDat$CHR38 <- paste0('chr', pidDat$CHR38)
  igadDat$CHR38 <- paste0('chr', igadDat$CHR38)

  pidGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  back_to_back_manhattan(pidGRanges, igadGRanges, top_label='PID', bottom_label='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-26e6', tick_dist=1e5, plot_genes=T)

  expect_equal(2, 1+1)
})

test_that("\'Zoomed\' Manhattan with three plots is drawn with genes", {
  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  sleDat <- data.table::fread(file=system.file('extdata', 'sleGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  pidDat <- subset(pidDat, -log10(P) < 15)
  igadDat <- subset(igadDat, -log10(P) < 15)
  sleDat <- subset(sleDat, -log10(P) < 15)

  pidDat$CHR38 <- paste0('chr', pidDat$CHR38)
  igadDat$CHR38 <- paste0('chr', igadDat$CHR38)
  sleDat$CHR38 <- paste0('chr', sleDat$CHR38)

  pidGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pidDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  igadGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igadDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
  sleGRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(sleDat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  back_to_back_manhattan(pidGRanges, igadGRanges, third_gRanges=sleGRanges, third_label='SLE', top_label='PID', bottom_label='IgAD', main='Chromosome 6, 25M-26M', zoom='chr6:25e6-26e6', tick_dist=1e5, plot_genes=T)

  expect_equal(2, 1+1)
})
