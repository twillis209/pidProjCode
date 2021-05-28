test_that("Error is thrown for missing p-value column with correct label when max_p is specified", {
  set.seed(42)
  gwas_dat <- data.table('CHR38' = 1:10, 'BP38' = 28567:28576, 'SNP' = paste0('rs', 1:10), 'not_p' = runif(10))

  expect_error(gwas_to_granges(gwas_dat, max_p = 0.1), "max_p specified but neither \'P\' nor \'p\' are names of metadata columns")
})

test_that("Error is thrown for missing p-value column with correct label when min_p is specified", {
  set.seed(42)
  gwas_dat <- data.table('CHR38' = 1:10, 'BP38' = 28567:28576, 'SNP' = paste0('rs', 1:10), 'not_p' = runif(10))

  expect_error(gwas_to_granges(gwas_dat, min_p = 0.1), "min_p specified but neither \'P\' nor \'p\' are names of metadata columns")
})

test_that("max_p works with p-value label \'P\'", {
  set.seed(42)
  gwas_dat <- data.table('CHR38' = 1:10, 'BP38' = 28567:28576, 'SNP' = paste0('rs', 1:10), 'P' = seq(0.1, 1, length.out = 10))

  expect_equal(mcols(gwas_to_granges(gwas_dat, max_p = 0.3))$P, c(0.1, 0.2))
})

test_that("max_p works with p-value label \'p\'", {
  set.seed(42)
  gwas_dat <- data.table('CHR38' = 1:10, 'BP38' = 28567:28576, 'SNP' = paste0('rs', 1:10), 'p' = seq(0.1, 1, length.out = 10))

  expect_equal(mcols(gwas_to_granges(gwas_dat, max_p = 0.3))$p, c(0.1, 0.2))
})

test_that("min_p works with p-value label \'P\'", {
  set.seed(42)
  gwas_dat <- data.table('CHR38' = 1:10, 'BP38' = 28567:28576, 'SNP' = paste0('rs', 1:10), 'P' = seq(0.1, 1, length.out = 10))

  expect_equal(mcols(gwas_to_granges(gwas_dat, min_p = 0.3))$P, c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
})

test_that("min_p works with p-value label \'p\'", {
  set.seed(42)
  gwas_dat <- data.table('CHR38' = 1:10, 'BP38' = 28567:28576, 'SNP' = paste0('rs', 1:10), 'p' = seq(0.1, 1, length.out = 10))

  expect_equal(mcols(gwas_to_granges(gwas_dat, min_p = 0.3))$p, c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
})

test_that("max_p and min_p work with p-value label \'P\'", {
  set.seed(42)
  gwas_dat <- data.table('CHR38' = 1:10, 'BP38' = 28567:28576, 'SNP' = paste0('rs', 1:10), 'P' = seq(0.1, 1, length.out = 10))

  expect_equal(mcols(gwas_to_granges(gwas_dat, min_p = 0.3, max_p = 0.7))$P, c(0.3, 0.4, 0.5, 0.6))
})

test_that("max_p and min_p work with p-value label \'p\'", {
  set.seed(42)
  gwas_dat <- data.table('CHR38' = 1:10, 'BP38' = 28567:28576, 'SNP' = paste0('rs', 1:10), 'p' = seq(0.1, 1, length.out = 10))

  expect_equal(mcols(gwas_to_granges(gwas_dat, min_p = 0.3, max_p = 0.7))$p, c(0.3, 0.4, 0.5, 0.6))
})
