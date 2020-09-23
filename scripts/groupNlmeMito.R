#' ---
#' title: "Fitting various NLME models for whole mitochondria"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(reshape2)
library(nlme)
library(plyr)
library(magrittr)

load("../intermediaryData/triexp.nlsList.LM.rda")

df_grouped <- read_experiments()

#' WT and Hacd1-KO together
#' ======================================================================================

nlme_outer_all <- fluoDecay_nlme(
  df_grouped,
  c("WT", "Hacd1-KO"), "mito",
  randomStr = c("B1", "B2", "B3", "theta1", "theta2", "theta3"),
  random = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1)
)
save(nlme_outer_all, file="../intermediaryData/mito.nlme_all.rda")

#' Allow both Bs and thetas to vary
#' ======================================================================================

pop_mitoplasts <- prime(triexp.nlsList, df_grouped, pop = c("WT", "Hacd1-KO"), type = "mito", c("B1", "B2", "B3", "theta1", "theta2", "theta3"))

randomStart <- pop_mitoplasts$randomStart
dat <- pop_mitoplasts$dat

meta <- ddply(df_grouped, c("File"), function (x) {
  data.frame(Population = x$Population[1], Type = x$Type[1])
})
meta_sub <- subset(meta, Type == "mito")
fixedStart <- c()
coefNames <- c("theta1", "B1", "theta2", "B2", "theta3", "B3")
for (coefName in coefNames) {
  a <- fixed.effects(nlme_outer_all)[[coefName]]
  b <- mean(random.effects(nlme_outer_all)[meta_sub$Population == "Hacd1-KO", coefName])
  c <- mean(random.effects(nlme_outer_all)[meta_sub$Population == "WT", coefName])
  fixedStart[paste0(coefName, ".(Intercept)")] <- a + b
  fixedStart[paste0(coefName, ".PopulationWT")] <- c - b
}

model <- triexpModelFactory(dat)
triexp_nlme_model <- formula(model)

nlme_outer_split <- nlme(
  triexp_nlme_model,
  fixed = theta1 + B1 + theta2 + B2 + theta3 + B3 ~ Population,
  random = list(
    File = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1)
  ),
  weights = varConstPower(fixed = c(power = 0.5), const = 1, form = ~ fitted(.)),
  data = dat %>% transform(Population = factor(Population, levels = c("Hacd1-KO", "WT"))),
  start = list(
    fixed = fixedStart,
    random = random.effects(nlme_outer_all)
  ),
  control = nlmeControl(
    pnlsTol = 1,
    msVerbose = TRUE
  ),
  verbose = TRUE
)
save(nlme_outer_split, file="../intermediaryData/mito.nlme_split.rda")

#' As a conclusion: nope, the two models are alike
print(anova(nlme_outer_all, nlme_outer_split))

#' Let's do the same with only Bs allowed to vary
#' ======================================================================================
pop_mitoplasts <- prime(triexp.nlsList, df_grouped, pop = c("WT", "Hacd1-KO"), type = "mito", c("B1", "B2", "B3", "theta1", "theta2", "theta3"))
dat <- pop_mitoplasts$dat

meta_sub <- subset(meta, Type == "mito")
fixedStart <- c()
fixedStart <- fixed.effects(nlme_outer_all)[c("theta1", "theta2", "theta3")]
coefNames <- c("B1", "B2", "B3")
for (coefName in coefNames) {
  a <- fixed.effects(nlme_outer_all)[[coefName]]
  b <- mean(random.effects(nlme_outer_all)[meta_sub$Population == "Hacd1-KO", coefName])
  c <- mean(random.effects(nlme_outer_all)[meta_sub$Population == "WT", coefName])
  fixedStart[paste0(coefName, ".(Intercept)")] <- a + b
  fixedStart[paste0(coefName, ".PopulationWT")] <- c - b
}

model <- triexpModelFactory(dat)
triexp_nlme_model <- formula(model)

nlme_outer_split_Bs <- nlme(
  triexp_nlme_model,
  fixed = list(theta1 + theta2 + theta3 ~ 1, B1 + B2 + B3 ~ Population),
  random = list(
    File = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1)
  ),
  weights = varConstPower(fixed = c(power = 0.5), const = 1, form = ~ fitted(.)),
  data = dat %>% transform(Population = factor(Population, levels = c("Hacd1-KO", "WT"))),
  start = list(
    fixed = fixedStart,
    random = random.effects(nlme_outer_all)
  ),
  control = nlmeControl(
    pnlsTol = 1,
    msVerbose = TRUE
  ),
  verbose = TRUE
)
save(nlme_outer_split_Bs, file="../intermediaryData/mito.nlme_split_Bs.rda")

#' Again, only splitting the Bs, we get the same result
print(anova(nlme_outer_all, nlme_outer_split_Bs, nlme_outer_split))

#' Now, let's split only the thetas
#' ======================================================================================
pop_mitoplasts <- prime(triexp.nlsList, df_grouped, pop = c("WT", "Hacd1-KO"), type = "mito", c("B1", "B2", "B3", "theta1", "theta2", "theta3"))
dat <- pop_mitoplasts$dat

meta_sub <- subset(meta, Type == "mito")
fixedStart <- c()
fixedStart <- fixed.effects(nlme_outer_all)[c("B1", "B2", "B3")]
coefNames <- c("theta1", "theta2", "theta3")
for (coefName in coefNames) {
  a <- fixed.effects(nlme_outer_all)[[coefName]]
  b <- mean(random.effects(nlme_outer_all)[meta_sub$Population == "Hacd1-KO", coefName])
  c <- mean(random.effects(nlme_outer_all)[meta_sub$Population == "WT", coefName])
  fixedStart[paste0(coefName, ".(Intercept)")] <- a + b
  fixedStart[paste0(coefName, ".PopulationWT")] <- c - b
}

model <- triexpModelFactory(dat)
triexp_nlme_model <- formula(model)

nlme_outer_split_thetas <- nlme(
  triexp_nlme_model,
  fixed = list(B1 + B2 + B3 ~ 1, theta1 + theta2 + theta3 ~ Population),
  random = list(
    File = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1)
  ),
  weights = varConstPower(fixed = c(power = 0.5), const = 1, form = ~ fitted(.)),
  data = dat %>% transform(Population = factor(Population, levels = c("Hacd1-KO", "WT"))),
  start = list(
    fixed = fixedStart,
    random = random.effects(nlme_outer_all)
  ),
  control = nlmeControl(
    pnlsTol = 1,
    msVerbose = TRUE
  ),
  verbose = TRUE
)
save(nlme_outer_split_thetas, file="../intermediaryData/mito.nlme_split_thetas.rda")

#' Interestingly, we do NOT get the same result here
print(anova(nlme_outer_all, nlme_outer_split_thetas, nlme_outer_split))
print(anova(nlme_outer_all, nlme_outer_split))

