#' ---
#' title: "NLS on individual experiments (with `nlsList` from package `nlme`)"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(plyr)
library(magrittr)
library(ggplot2)
library(nlme)

df_grouped <- read_experiments()

df_grouped_all <- df_grouped
meta <- ddply(df_grouped_all, c("File"), function (x) {
  data.frame(Population = x$Population[1], Type = x$Type[1])
})
# This comes from a prior NLS fit on all the data together (data not shown)
start <- structure(c(1.79953624594729, 0.0753975243025479, -1.05834556059257,
  0.00875804362445855, -0.243275701393101, 0.0184522361929007), .Names = c("theta1",
    "B1", "theta2", "B2", "theta3", "B3"))

model <- triexpModelFactory(df_grouped_all)
triexp_nlme_model <- formula(model)

st <- proc.time()
triexp.nlsList2 <- nlsList(
  triexp_nlme_model,
  data = df_grouped_all,
  start = start
)
print(proc.time() - st)

# Let's save this
save(triexp.nlsList2, file="../intermediaryData/triexp.nlsList.rda")

#' We can plot intervals
plot(intervals(triexp.nlsList2))

#' You can see that for a number of experiments, convergence did NOT happen:
no_conv <- apply(intervals(triexp.nlsList2), 1, function (x) any(is.na(x)))
print(names(which(no_conv)))
print(sum(no_conv))

#' Nota: by using `nlsLM` with weights, this goes a little bit better
#' (see simpleIntervals2.R).
#'
#' Thetas are not necessarily ordered from low to high (i.e. theta1 lowest),
#' Which is problematic for comparison
thetas_order <- t(apply(
  coef(triexp.nlsList2),
  1,
  function (x) {
    thetas <- x[c("theta1", "theta2", "theta3")]
    setNames(order(thetas), c("low", "mid", "high"))
  }
))
print(thetas_order)

#' We can get a feeling of the result by plotting on a tau/B1 axis
#' To get tau from theta, with confidence intervals, we need to do a
#' transformation with some Monte Carlo sampling
triexp.vcovsMeans <- vcovsMeans(triexp.nlsList2)
df_freq_plot2 <- freqPlotCov(triexp.vcovsMeans)
print(plotFreq(df_freq_plot2))

#' We see that the higher the tau, the more difficult it is to get
#' to a precise pinning of it
#'
#' Now, let's plot the means and the variances of the means for
#' the interaction of population (WT/Hacd1-KO) and type (mito/mitoplaste).
#' To do this we must re-order the thetas and Bs.
triexp.varsMeans.grouped <- group.varsMeans(meta, triexp.vcovsMeans)

df_freq_plot_grouped <- freqPlotVar(triexp.varsMeans.grouped)
print(plotFreq(df_freq_plot_grouped))

#' In the end, we see a lack of significativity
#' that could be addressed with
#' 1) weighing and Levenberg–Marquardt fitting (simpleIntervals2);
#' 2) mixed effects modelling (groupNlme).

