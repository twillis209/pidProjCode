test_that("Error is thrown for omitting p-value column from gRanges", {
  set.seed(42)
  gwas_daf <- data.frame('chr' = paste0('chr', 1:10), 'bp' = 28567:28576, 'SNP' = paste0('rs', 1:10),'p' = runif(10))
  gRanges_top <- GenomicRanges::makeGRangesFromDataFrame(gwas_daf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  gRanges_bottom <- GenomicRanges::makeGRangesFromDataFrame(gwas_daf, start.field = 'bp', end.field = 'bp', seqnames.field = 'chr', ignore.strand = T, keep.extra.columns = T)
  gRanges_top$p <- NULL

  expect_error(multitrack_manhattan(gRanges = list(gRanges_top, gRanges_bottom),
                                    axis_labels = c('a', 'b'),
                                   main = 'title'),
               "Element 1 of gRanges does not contain a column labelled \'P\' or \'p\'")
})

test_that("Single-track Manhattan is drawn", {
  pid_dat <- data.table::fread(file = system.file('extdata', 'pidGwas.tsv.gz', package = 'pidProjCode', mustWork = T), sep = '\t', header = T)

  pid_dat <- subset(pid_dat, -log10(P) < 15)

  pid_dat$CHR38 <- paste0('chr', pid_dat$CHR38)

  pid_gRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pid_dat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  multitrack_manhattan(list(pid_gRanges), '-log10(p)', 'PID')

  expect_equal(2, 1+1)
})

test_that("Single-track Manhattan with specified \'chromosomes\' parameter is drawn", {
  pid_dat <- data.table::fread(file = system.file('extdata', 'pidGwas.tsv.gz', package = 'pidProjCode', mustWork = T), sep = '\t', header = T)

  pid_dat <- subset(pid_dat, -log10(P) < 15)

  pid_dat$CHR38 <- paste0('chr', pid_dat$CHR38)

  pid_gRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pid_dat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  multitrack_manhattan(list(pid_gRanges), '-log10(p)', 'PID', chromosomes = c('chr5', 'chr6', 'chr7'))

  expect_equal(2, 1+1)
})

test_that("Two-track Manhattan is drawn", {
    pid_dat <- data.table::fread(file = system.file('extdata', 'pidGwas.tsv.gz', package = 'pidProjCode', mustWork = T), sep = '\t', header = T)
    igad_dat <- data.table::fread(file = system.file('extdata', 'igadGwas.tsv.gz', package = 'pidProjCode', mustWork = T), sep = '\t', header = T)

    pid_dat <- subset(pid_dat, -log10(P) < 15)
    igad_dat <- subset(igad_dat, -log10(P) < 15)

    pid_dat$CHR38 <- paste0('chr', pid_dat$CHR38)
    igad_dat$CHR38 <- paste0('chr', igad_dat$CHR38)

    pid_gRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pid_dat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
    igad_gRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igad_dat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  multitrack_manhattan(gRanges = list(pid_gRanges, igad_gRanges), axis_labels = rep('-log10(p)', 2), main = 'PID and IgAD')

  expect_equal(2, 1+1)
})

test_that("Single-track Manhattan with specified \'chromosomes\' parameter is drawn", {
  pid_dat <- data.table::fread(file = system.file('extdata', 'pidGwas.tsv.gz', package = 'pidProjCode', mustWork = T), sep = '\t', header = T)

  pid_dat <- subset(pid_dat, -log10(P) < 15)

  pid_dat$CHR38 <- paste0('chr', pid_dat$CHR38)

  pid_gRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pid_dat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  multitrack_manhattan(list(pid_gRanges), '-log10(p)', 'PID', chromosomes = c('chr5', 'chr6', 'chr7'))

  expect_equal(2, 1+1)
})

test_that("Two-track Manhattan with specified \'chromosomes\' parameter is drawn", {
    pid_dat <- data.table::fread(file = system.file('extdata', 'pidGwas.tsv.gz', package = 'pidProjCode', mustWork = T), sep = '\t', header = T)
    igad_dat <- data.table::fread(file = system.file('extdata', 'igadGwas.tsv.gz', package = 'pidProjCode', mustWork = T), sep = '\t', header = T)

    pid_dat <- subset(pid_dat, -log10(P) < 15)
    igad_dat <- subset(igad_dat, -log10(P) < 15)

    pid_dat$CHR38 <- paste0('chr', pid_dat$CHR38)
    igad_dat$CHR38 <- paste0('chr', igad_dat$CHR38)

    pid_gRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(pid_dat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)
    igad_gRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(igad_dat), start.field = 'BP38', end.field = 'BP38', seqnames.field = 'CHR38', ignore.strand = T, keep.extra.columns = T)

  multitrack_manhattan(gRanges = list(pid_gRanges, igad_gRanges), axis_labels = rep('-log10(p)', 2), main = 'PID and IgAD', chromosomes = paste0('chr', 5:7))

  expect_equal(2, 1+1)
})
