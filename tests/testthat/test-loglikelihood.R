library(survtmle)
context("Testing computation of the log-likelihood of a given model.")

test_that("log-likelihood fails when invoked improperly", {

  expect_error( LogLikelihood(Beta = ,
                              X = ,
                              Y = )
              )

})
