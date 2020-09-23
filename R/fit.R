#' Perform an `nlme` or `gnls` depending on the parameters given.
#'
#' @param df_grouped Dataset
#' @param pop To subset data to some population levels
#' @param type To subset data to some type levels
#' @param randomStr Character vector of parameter names that have random effects
#' @param random Random effect formula (see `nlme`)
#' @param fixed Fixed effect formula (see `nlme`)
#' @param weights Whether to have const+power weights
#' @param control `nlmeControl` object
#'
#' @note
#' If `randomStr == NULL`, a `gnls` is performed.
#'
#' We pollute the global environment because of the way R scopes work for formulas.
fluoDecay_nlme <- function (df_grouped, pop, type, randomStr, random, fixed = theta1 + B1 + theta2 + B2 + theta3 + B3 ~ 1, weights = TRUE, control = nlmeControl(msVerbose = TRUE)) {
  # NOTE: We need globals: that's why R's scope is actually horrible.
  pop_WT_mito <<- prime(triexp.nlsList, df_grouped, pop, type, randomStr)

  fixedStart <<- pop_WT_mito$fixedStart
  randomStart <<- pop_WT_mito$randomStart
  dat <<- pop_WT_mito$dat

  model <<- triexpModelFactory(dat)
  triexp_nlme_model <<- formula(model)

  start <- if (is.null(randomStr)) {
    fixedStart
  } else {
    list(
      fixed = fixedStart,
      random = randomStart
    )
  }

  st <- proc.time()
  args <- list(
    triexp_nlme_model,
    fixed = fixed,
    random = random,
    data = dat,
    start = start,
    weights = if (!weights) {
      NULL
    } else {
      varConstPower(fixed = c(power = 0.5), const = 1, form = ~ fitted(.))
    },
    verbose = TRUE,
    control = control
  )
  if (is.null(random)) {
    args$random <- NULL
    args$fixed <- NULL
    fit <- do.call(gnls, args)
  } else {
    fit <- do.call(nlme, args)
  }
  print(proc.time() - st)

  fit
}
