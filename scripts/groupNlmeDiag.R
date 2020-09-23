#' ---
#' title: "Determining the \"class\" of NLME models to use"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

#' Follow-up from `simpleIntervals.R`.
#'
#' The first thing we do is to figure out what NLME model we need to use.
#' We tried different configurations and tested them using ANOVA.
#'
#' We ended up using a diagonal random effects configuration with
#' `const + power` intra-experiment variance structure.
#'
#' The `const + power` variance structure is a great way to model Poisson
#' distributions (e.g., distribution of event counts) without going into
#' generalized equations.
#'
#' The diagonal random-effect structure is somewhat a default.  Using
#' more or less refined models necessitate a justification, which are given below.

library(Prola2020.Sci.Adv)
library(nlme)

load("../intermediaryData/triexp.nlsList.LM.rda")

df_grouped <- read_experiments()

#' Without weights nor random effect: do not converge
#' ======================================================================================
withErr <- FALSE
tryCatch({
  invisible(fluoDecay_nlme(df_grouped, "WT", "mitoplaste", randomStr = NULL, random = NULL, weights = FALSE))
}, error = function (e) {
  print(e$message)
  withErr <<- TRUE
})
if (!withErr) { stop("error expected") }

#' With weights only: do not converge
#' ======================================================================================
withErr <- FALSE
tryCatch({
  invisible(fluoDecay_nlme(df_grouped, "WT", "mito", randomStr = NULL, random = NULL, weights = TRUE))
}, error = function (e) {
  print(e$message)
  withErr <<- TRUE
})
if (!withErr) { stop("error expected") }

#' Let's see with full diagonal
#' ======================================================================================
nls_mito_diag <- fluoDecay_nlme(df_grouped, "WT", "mitoplaste", randomStr = c("B1", "B2", "B3", "theta1", "theta2", "theta3"), random = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1))

#' Let's see with half diagonals
#' ======================================================================================
nls_mito_diag_Bs <- fluoDecay_nlme(df_grouped, "WT", "mitoplaste", randomStr = c("B1", "B2", "B3"), random = pdDiag(B1 + B2 + B3 ~ 1))
nls_mito_diag_thetas <- fluoDecay_nlme(df_grouped, "WT", "mitoplaste", randomStr = c("theta1", "theta2", "theta3"), random = pdDiag(theta1 + theta2 + theta3 ~ 1))

#' Full diagonal and off-diagonal terms only for B1/B2/B3 (because they resolve to proportions)
#' ======================================================================================
withErr <- FALSE
tryCatch({
  fluoDecay_nlme(
    df_grouped,
    "WT", "mitoplaste",
    randomStr = c("B1", "B2", "B3", "theta1", "theta2", "theta3"),
    random = pdBlocked(list(B1 + B2 + B3 ~ 1, pdDiag(theta1 + theta2 + theta3 ~ 1)))
  )
}, error = function (e) {
  print(e$message)
  withErr <<- TRUE
})
if (!withErr) { stop("error expected") }

#' Conclusion
#' ======================================================================================

#' Anova full/half diagonals vs. weights only.
#' The full diagonal model is clearly better.
#' Notice that what brings value is mainly variability in thetas, not Bs.
#' But since we cannot find a biological justification for this, we stick
#' to the full diagonal model.
print(anova(nls_mito_diag, nls_mito_diag_Bs, nls_mito_diag_thetas))

save(
  nls_mito_diag,
  nls_mito_diag_Bs,
  nls_mito_diag_thetas,
  file = "../intermediaryData/groupNlmeDiag.rda"
)

