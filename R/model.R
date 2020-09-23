#' Generate a model function from a dataset (you need the dataset to get the prompt).
#'
#' @param dat A `data.frame` with the following columns (and possibly other, ignored
#' columns): `File` (a factor to designate the source file/experiment), `Chan`
#' (the measurement channels, for instance from 1 to 4096), `Prompt` (the instrumental
#' response, or "prompt", for the corresponding measurement channel), `Decay` (the
#' experimental response, or "decay", for the corresponding measurement channel).
#' @param calcGradient Whether to calculate the gradient and pass it as an attribute
#' to the return value (`TRUE`, default) or to skip this step (`FALSE`).
#' @param multithread Whether to do the computation on multiple threads (on the C++ side,
#' R is monothreaded).
#' @return A function that returns a numeric vector of calculated experimental
#' responses from the following arguments: `Chan` (the channels for which to calculate
#' the response), `File` (the files for which ...), `thetas` parameters (resolving to
#' characteristic times), `Bs` parameters (scaling parameters), `timeinterval` (a
#' numeric value giving the time interval between two channels; this allows us to give
#' proper scaling for `thetas`).
#'
#' @notes
#' If you need a time shift or to add a nonzero background noise level, you need to
#' implement that on your own.  This can be done directly in the NLS formula.
#'
#' The `thetas` and `Bs` parameters are either one-rowed matrices or matrices with the same
#' number of rows as there are rows in the dataset.  The number of columns give the number
#' of exponentials.
#'
#' More info on the model and meaning of the parameters can be found in
#' `scripts/article/supplementaryMaterials.pdf`.
baseModelFactory <- function (dat, calcGradient = TRUE, multithread = TRUE) {
  # Get prompt starts, ends and the maximum Chan for each experiment
  dat$Row <- 1:NROW(dat)
  promptStarts <- sapply(levels(dat$File), function (file) min(subset(dat, File == file)$Row))
  promptEnds <- sapply(levels(dat$File), function (file) max(subset(dat, File == file)$Row))
  iiMax <- sapply(levels(dat$File), function (file) {
    ss <- subset(dat, File == file)
    iMax <- max(ss$Chan)
    if (iMax != NROW(ss) || any(duplicated(ss$Chan))) {
      stop("Channels must be contiguous, start from 1 and not be duplicated")
    }
    iMax
  })

  # Get parameters actually given to the C++ implementation from the
  # parameters given to the R model function.
  modelParams <- function (Chan, File, thetas, Bs, timeinterval) {
    # Validate some of the assumptions
    if (NCOL(thetas) != NCOL(Bs)) stop("thetas and Bs must be of same length")
    if (NROW(thetas) != 1 && NROW(thetas) != length(Chan)) {
      stop("'thetas' should either be of length 1 or the number of channels")
    }
    if (NROW(Bs) != 1 && NROW(Bs) != length(Chan)) {
      stop("'Bs' should either be of length 1 or the number of channels")
    }

    if (NROW(thetas) != 1) {
      # If thetas and Bs are of length the number of channels, we gather unique sets of parameters.
      # This allows us to call the underlying C++ objective function less (we do this once
      # per unit set per file) and more importantly to do a Fast Fourier Transform.
      # We also bind `thetas` and `Bs` into a unique `params` matrix.
      params <- cbind(thetas, Bs, as.numeric(File))

      # To get to unique parameters, we do not use unique, but a fast C++ implementation
      # that is around a hundred times faster (this particular line of code rapidly became
      # a bottleneck).
      uniqueParams <- fastNumericUnique(params, 10)

      # BTW, the reason we need to sort out unique parameters is due to the way NLME
      # calls the model.
    } else {
      # Otherwise, we reuse the only set of parameters across the files.
      nlvl <- nlevels(File)
      params <- cbind(
        matrix(
          ncol = NCOL(thetas) * 2,
          byrow = TRUE,
          rep(c(thetas, Bs), nlvl)
        ),
        1:nlvl
      )
      uniqueParams <- params
    }

    # We loop over each unique set of parameters and files and create plans.
    # Those plans are executed in parallel in C++ (truly, with threads, no forking,
    # something impossible to do in R which is not multithreaded)
    #
    # NOTE: We could probably do this directly in the C++ implementation, but
    # the impact on performance is really negligible.
    plans <- list()
    for (m in 1:NROW(uniqueParams)) {
      file = levels(File)[uniqueParams[m, NCOL(uniqueParams)]]
      rowSubset <-
        if (NROW(thetas) == 1) {
          TRUE
        } else {
          colSums(t(params) == uniqueParams[m,]) == NCOL(uniqueParams)
        }
      rowSubset2 <- rowSubset & File == file
      whichRowSubset <- which(rowSubset2)
      plans[[length(plans) + 1]] <- list(
        rowSubset = whichRowSubset - 1,
        m = m,
        file = as.numeric(factor(file, levels(dat$File))),
        iiMax = iiMax[[file]]
      )
    }

    if (length(plans) == 0) {
      stop("no plan!")
    }

    list(
      length(Chan),
      Chan,
      t(uniqueParams[, 1:NCOL(thetas), drop=F]),
      t(uniqueParams[, (NCOL(thetas) + 1):(NCOL(thetas) * 2), drop=F]),
      timeinterval = timeinterval,
      dat$Prompt,
      promptStarts,
      promptEnds,
      plans,
      calcGradient,
      multithread
    )
  }

  model <- function (...) {
    params <- modelParams(...)
    do.call(decayModelExecute, params)
  }

  return (model)
}
