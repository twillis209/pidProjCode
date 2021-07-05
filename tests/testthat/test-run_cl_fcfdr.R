test_that("Error is thrown when aux_transform length does not match length of auxiliary", {
  cl_args <- c("-p", "P.pid", "-q", "P.ra,P.aster", "-w", "Weight.ra", "-v", "v.ra", "-op", "test.tsv.gz", "-nt", "8", "-at", "log")

  expect_error(run_cl_fcfdr(cl_args), "No. of specified auxiliary transformations does not match number of auxiliary covariates")
})

test_that("Error is thrown when aux_transform value is invalid", {
  cl_args <- c("-p", "P.pid", "-q", "P.ra", "-w", "Weight.ra", "-v", "v.ra", "-op", "test.tsv.gz", "-nt", "8", "-at", "invalid_value")

  expect_error(run_cl_fcfdr(cl_args), "Specified auxiliary transformations contain values not matching one of following valid transformation identifiers: \'identity\', \'log\', \'z\'")
})
