library(testthat)

test_that("The surprisal analysis function works for a valid csv file.", {

  data <- read.csv(system.file("extdata", "sample_expression.csv", package = "SurprisalAnalysis"), header=TRUE)

  results <- surprisal_analysis(data)

  expect_length(results, 2)

  expect_type(results[[1]], "double")

  expect_length(results[[1]], (ncol(data)-1)^2)

  expect_type(results[[2]], "double")

  expect_length(results[[2]], (ncol(data)-1)*nrow(data))

})
