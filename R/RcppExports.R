# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' The `Worker` paradigm comes from `RcppParallel`.
NULL

#' Fast extraction of unique rows of a numeric matrix
#' @param expected: expected number of unique rows to pre-reserve memory and avoid costly
#' resizing.
fastNumericUnique <- function(x, expected) {
    .Call('_Prola2020_Sci_Adv_fastNumericUnique', PACKAGE = 'Prola2020.Sci.Adv', x, expected)
}

decayModelExecute <- function(outputLength, ii, allAlphas, allBetas, timeinterval, prompts, promptStarts, promptEnds, plans, grad, multithread) {
    .Call('_Prola2020_Sci_Adv_decayModelExecute', PACKAGE = 'Prola2020.Sci.Adv', outputLength, ii, allAlphas, allBetas, timeinterval, prompts, promptStarts, promptEnds, plans, grad, multithread)
}

