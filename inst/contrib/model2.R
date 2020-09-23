library(reshape2)
library(plyr)

# Let's have some generic data
timeinterval <- 0.02743484
set.seed(123)

# We set the prompt to be a gamma distribution
prompt <- function () {
  dgamma(1:4096, shape = runif(1, min=3.8, max=4.2), scale = runif(1, min=90, max=110)) * 100
}

# Parameters of our models
parameters <- function (nExponentials) {
  list(
    thetas = sapply(
      -seq(from = -0.25, by = 1, length.out = nExponentials),
      function (mu) rnorm(1, mu, 0.05)
    ),
    Bs = runif(nExponentials, 0.0001, 0.05)
  )
}

calcDecay <- function (Chan, Prompt, pars) {
  # Get non-convolved version
  beforeConv <- vector(mode = "numeric", length = length(Chan))
  for (i in 1:length(pars$thetas)) {
    beforeConv <- beforeConv +
      pars$Bs[[i]] * exp(-exp(pars$thetas[[i]]) * Chan * timeinterval)
  }

  # Convolve.  Voluntarily, we do not use `convolve` from `utils` (it relies
  # on a FFT).
  sapply(
    Chan,
    function (i) {
      sum(beforeConv[1:i] * Prompt[i - 1:i + 1])
    }
  )
}

fastSlowModelTest <- function (nExponentials) {
  # Files
  files <- paste0("file_", 1:5)
  params <- setNames(
    replicate(length(files), { parameters(nExponentials) }, simplify = FALSE),
    files
  )

  df_test_dataset <- rbind.fill(lapply(
    files,
    function (file) {
      data.frame(File = file, Chan = 1:4096, Prompt = prompt()) %>%
        transform(Decay = calcDecay(Chan, Prompt, params[[file]]))
    }
  ))
  df_test_dataset$File <- factor(df_test_dataset$File, files)
  df_test_dataset$Row <- 1:NROW(df_test_dataset)
  df_test_dataset$Time = (df_test_dataset$Chan - 1) * timeinterval

  model <- baseModelFactory(df_test_dataset)
  fastDecay <- with(df_test_dataset, model(
    Chan, File,
    t(sapply(params[File], function (x) x$thetas)),
    t(sapply(params[File], function (x) x$Bs)),
    timeinterval
  ))

  # Compute the numerical Jacobian to compare it to the calculated Jacobian
  numGrad <- do.call(rbind, lapply(files, function (file) {
    pars <- params[[file]]
    around <- c()
    for (i in 1:length(pars$thetas)) {
      around[paste0("theta", i)] <- pars$thetas[[i]]
      around[paste0("B", i)] <- pars$Bs[[i]]
    }
    dat <- subset(df_test_dataset, File == file)
    numDeriv::jacobian(
      function (x) {
        thetas <- x[grepl("^theta", names(x))]
        Bs <- x[grepl("^B", names(x))]
        calcDecay(dat$Chan, dat$Prompt, list(thetas = thetas, Bs = Bs))
      },
      around,
      method = "simple"
    )
  }))
  colnames(numGrad) <- Reduce(function (acc, n) c(acc, paste0("theta", n), paste0("B", n)), 1:nExponentials, c())
  calcGrad <- attr(fastDecay, "gradient")

  df_gradient <- merge(
    melt(numGrad, varnames=c("Row", "variable"), value.name="Num"),
    melt(calcGrad, varnames=c("Row", "variable"), value.name="Calc")
  )

  list(
    df_dataset = df_test_dataset %>%
      transform(FastDecay = as.numeric(fastDecay)),
    df_gradient = df_gradient
  )
}
