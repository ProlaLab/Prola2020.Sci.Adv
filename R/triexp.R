#' Return a triexponential model from a dataset.
triexpModelFactory <- function (dat) {
  baseModel <- baseModelFactory(dat)

  structure(
    function (Chan, File, theta1, theta2, theta3, B1, B2, B3) {
      baseModel(Chan, File, cbind(theta1, theta2, theta3), cbind(B1, B2, B3), timeinterval)
    },
    class = c("function", "triexpModel")
  )
}

#' Get the formula of the triexponential model, to be plugged in `nls`, `nlme`, ...
formula.triexpModel <- function (model) {
  modelFnName <- match.call()$model
  formula(paste0(
    "Decay ~ ",
    modelFnName,
    "(Chan, File, theta1, theta2, theta3, B1, B2, B3)"
  ), env = .GlobalEnv)
}
