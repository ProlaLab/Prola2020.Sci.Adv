context("model of fluorescence decay")

source(system.file("contrib/model2.R", package="Prola2020.Sci.Adv"), chdir=T)

testFastSlowModels <- function (ans) {
  df_dataset <- ans$df_dataset
  df_gradient <- ans$df_gradient

  # We expect fast and slow decay to differ by less than 3% after the initial divergence
  # (due to numerical instability).
  expect_lt(max(abs((df_dataset$Decay - df_dataset$FastDecay) / df_dataset$Decay)[df_dataset$Time > 5]), 0.03, 'Decay difference')

  # We expect fast and slow decay to have a less than 1% difference in the gradient
  expect_lt(mean(abs(df_gradient$Num - df_gradient$Calc)) / mean(abs(df_gradient$Num)), 0.01, 'Gradient difference')
}

test_that("the test dataset and the C++ model agree in the base case", {
  options(warn = 2)
  ans <- fastSlowModelTest(nExponentials = 3)
  testFastSlowModels(ans)
})

test_that("when there are four exponentials", {
  options(warn = 2)
  ans <- fastSlowModelTest(nExponentials = 4)
  testFastSlowModels(ans)
})
