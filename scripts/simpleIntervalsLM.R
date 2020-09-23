#' ---
#' title: "NLS on individual experiments (with `nlsLM` from `minpack.lm` - Levenberg-Marquardt NLS)"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)

library(ggplot2)
library(parallel)
library(plyr)
library(magrittr)
library(minpack.lm)

df_grouped <- read_experiments()

df_grouped_all <- df_grouped
files <- unique(df_grouped_all$File)
meta <- ddply(df_grouped_all, c("File"), function (x) {
  data.frame(Population = x$Population[1], Type = x$Type[1])
})
start <- structure(c(1.79953624594729, 0.0753975243025479, -1.05834556059257,
  0.00875804362445855, -0.243275701393101, 0.0184522361929007), .Names = c("theta1",
    "B1", "theta2", "B2", "theta3", "B3"))

triexp.nlsList <- list()
st <- proc.time()
for (file in files) {
  dat <- subset(df_grouped, File == file)
  model <- triexpModelFactory(dat)
  triexp_nlme_model <- formula(model)
  cat("***", file, "***\n")
  ans <- nlsLM(
    triexp_nlme_model,
    data = dat,
    start = start,
    weights = sqrt(dat$Decay),
    lower = c(-Inf, 1e-5, -Inf, 1e-5, -Inf, 1e-5)
  )

  triexp.nlsList[[file]] <- ans
}
print(proc.time() - st)

#' Let's save this
save(triexp.nlsList, file="../intermediaryData/triexp.nlsList.LM.rda")

#' Interestingly, with LM, we get convergence every time
#' The intervals are as follows:
triexp.vcovsMeans <- vcovsMeans(triexp.nlsList)
df_freq_plot2 <- freqPlotCov(triexp.vcovsMeans)
print(plotFreq(df_freq_plot2))

#' We see here that there is a blow-up for one of the Hacd1-KO.mitoplaste experiment.
#' Let's censor this blow-up.
censoredFiles <- c("mitoplaste 8 C filter")
triexp.vcovsMeans.censored <- triexp.vcovsMeans
triexp.vcovsMeans.censored[[censoredFiles]] <- NULL
df_freq_plot2 <- freqPlotCov(triexp.vcovsMeans.censored)
print(plotFreq(df_freq_plot2))

#' Let's group by population/type
triexp.varsMeans.grouped <- group.varsMeans(subset(meta, !(File %in% censoredFiles)), triexp.vcovsMeans.censored)
df_freq_plot_grouped <- freqPlotVar(triexp.varsMeans.grouped)
print(plotFreq(df_freq_plot_grouped))
