test_that("Sample is drawn without error in simple use case", {
  testDaf <- data.frame(a=1:6, b=7:12)
  row.names(testDaf) <- c('a', 'b', 'c', 'd', 'e', 'f')
  expect_equal(prettyPrintDataFrame(testDaf),
               '| |a|b|\n|a|1|7|\n|b|2|8|\n|c|3|9|\n|d|4|10|\n|e|5|11|\n|f|6|12|\n')
})
