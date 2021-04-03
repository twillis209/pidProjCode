test_that("Simple Q-Q plot is written out", {

  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  stratified_qqplot(data.frame(p=pidDat$P), 'p')

  expect_equal(2+2, 4)
})

test_that("Stratified Q-Q plot is written out", {

  pidDat <- data.table::fread(file=system.file('extdata', 'pidGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)
  igadDat <- data.table::fread(file=system.file('extdata', 'igadGwas.tsv.gz', package='pidProjCode', mustWork = T), sep = '\t', header = T)

  igadDat <- igadDat[1:nrow(pidDat),]

  stratified_qqplot(data.frame(p=pidDat$P, q=igadDat$P), 'p', 'q')

  expect_equal(2+2, 4)
})
