#' ---
#' title: "Fitting various NLME models for mitoplasts"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(reshape2)
library(nlme)
library(magrittr)
library(plyr)

load("../intermediaryData/triexp.nlsList.LM.rda")

df_grouped <- read_experiments()

#' WT and Hacd1-KO together
#' ======================================================================================
#'
#' This censorship comes from analysis in `simpleIntervals2.R`.
censoredFiles <- c("mitoplaste 8 C filter")
nlme_all <- fluoDecay_nlme(
  df_grouped %>% subset(!(File %in% censoredFiles)),
  c("WT", "Hacd1-KO"), "mitoplaste",
  randomStr = c("B1", "B2", "B3", "theta1", "theta2", "theta3"),
  random = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1)
)
save(nlme_all, file="../intermediaryData/mitoplasts.nlme_all.rda")

#' Allow both Bs and thetas to vary
#' ======================================================================================

pop_mitoplasts <- prime(triexp.nlsList, df_grouped %>% subset(!(File %in% censoredFiles)), pop = c("WT", "Hacd1-KO"), type = "mitoplaste", c("B1", "B2", "B3", "theta1", "theta2", "theta3"))

randomStart <- pop_mitoplasts$randomStart
dat <- pop_mitoplasts$dat

meta <- ddply(df_grouped, c("File"), function (x) {
  data.frame(Population = x$Population[1], Type = x$Type[1])
})
meta_sub <- subset(meta, Type == "mitoplaste" & !(File %in% censoredFiles))
fixedStart <- c()
coefNames <- c("theta1", "B1", "theta2", "B2", "theta3", "B3")
for (coefName in coefNames) {
  a <- fixed.effects(nlme_all)[[coefName]]
  b <- mean(random.effects(nlme_all)[meta_sub$Population == "Hacd1-KO", coefName])
  c <- mean(random.effects(nlme_all)[meta_sub$Population == "WT", coefName])
  fixedStart[paste0(coefName, ".(Intercept)")] <- a + b
  fixedStart[paste0(coefName, ".PopulationWT")] <- c - b
}

model <- triexpModelFactory(dat)
triexp_nlme_model <- formula(model)

nlme_split <- nlme(
  triexp_nlme_model,
  fixed = theta1 + B1 + theta2 + B2 + theta3 + B3 ~ Population,
  random = list(
    File = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1)
  ),
  weights = varConstPower(fixed = c(power = 0.5), const = 1, form = ~ fitted(.)),
  data = dat %>% transform(Population = factor(Population, levels = c("Hacd1-KO", "WT"))),
  start = list(
    fixed = fixedStart,
    random = random.effects(nlme_all)
  ),
  control = nlmeControl(
    pnlsTol = 1,
    msVerbose = TRUE
  ),
  verbose = TRUE
)
save(nlme_split, file="../intermediaryData/mitoplasts.nlme_split.rda")

#' As a conclusion: yes, there is something, the two models are different
print(anova(nlme_all, nlme_split))

#' Now, let's try to do the same with only the Bs to vary
#' ======================================================================================
pop_mitoplasts <- prime(triexp.nlsList, df_grouped %>% subset(!(File %in% censoredFiles)), pop = c("WT", "Hacd1-KO"), type = "mitoplaste", c("B1", "B2", "B3", "theta1", "theta2", "theta3"))
dat <- pop_mitoplasts$dat

meta_sub <- subset(meta, Type == "mitoplaste" & !(File %in% censoredFiles))
fixedStart <- c()
fixedStart <- fixed.effects(nlme_all)[c("theta1", "theta2", "theta3")]
coefNames <- c("B1", "B2", "B3")
for (coefName in coefNames) {
  a <- fixed.effects(nlme_all)[[coefName]]
  b <- mean(random.effects(nlme_all)[meta_sub$Population == "Hacd1-KO", coefName])
  c <- mean(random.effects(nlme_all)[meta_sub$Population == "WT", coefName])
  fixedStart[paste0(coefName, ".(Intercept)")] <- a + b
  fixedStart[paste0(coefName, ".PopulationWT")] <- c - b
}

model <- triexpModelFactory(dat)
triexp_nlme_model <- formula(model)

nlme_split_Bs <- nlme(
  triexp_nlme_model,
  fixed = list(theta1 + theta2 + theta3 ~ 1, B1 + B2 + B3 ~ Population),
  random = list(
    File = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1)
  ),
  weights = varConstPower(fixed = c(power = 0.5), const = 1, form = ~ fitted(.)),
  data = dat %>% transform(Population = factor(Population, levels = c("Hacd1-KO", "WT"))),
  start = list(
    fixed = fixedStart,
    random = random.effects(nlme_all)
  ),
  control = nlmeControl(
    pnlsTol = 1,
    msVerbose = TRUE
  ),
  verbose = TRUE
)
save(nlme_split_Bs, file="../intermediaryData/mitoplasts.nlme_split_Bs.rda")

#' Again, only splitting the Bs, we get the same result
print(anova(nlme_all, nlme_split_Bs, nlme_split))

#' Now, let's split only the thetas
#' ======================================================================================
pop_mitoplasts <- prime(triexp.nlsList, df_grouped %>% subset(!(File %in% censoredFiles)), pop = c("WT", "Hacd1-KO"), type = "mitoplaste", c("B1", "B2", "B3", "theta1", "theta2", "theta3"))
dat <- pop_mitoplasts$dat

meta_sub <- subset(meta, Type == "mitoplaste" & !(File %in% censoredFiles))
fixedStart <- c()
fixedStart <- fixed.effects(nlme_all)[c("B1", "B2", "B3")]
coefNames <- c("theta1", "theta2", "theta3")
for (coefName in coefNames) {
  a <- fixed.effects(nlme_all)[[coefName]]
  b <- mean(random.effects(nlme_all)[meta_sub$Population == "Hacd1-KO", coefName])
  c <- mean(random.effects(nlme_all)[meta_sub$Population == "WT", coefName])
  fixedStart[paste0(coefName, ".(Intercept)")] <- a + b
  fixedStart[paste0(coefName, ".PopulationWT")] <- c - b
}

model <- triexpModelFactory(dat)
triexp_nlme_model <- formula(model)

nlme_split_thetas <- nlme(
  triexp_nlme_model,
  fixed = list(B1 + B2 + B3 ~ 1, theta1 + theta2 + theta3 ~ Population),
  random = list(
    File = pdDiag(B1 + B2 + B3 + theta1 + theta2 + theta3 ~ 1)
  ),
  weights = varConstPower(fixed = c(power = 0.5), const = 1, form = ~ fitted(.)),
  data = dat %>% transform(Population = factor(Population, levels = c("Hacd1-KO", "WT"))),
  start = list(
    fixed = fixedStart,
    random = random.effects(nlme_all)
  ),
  control = nlmeControl(
    pnlsTol = 1,
    msVerbose = TRUE
  ),
  verbose = TRUE
)
save(nlme_split_thetas, file="../intermediaryData/mitoplasts.nlme_split_thetas.rda")

#' Interestingly, we do NOT get the same result here
print(anova(nlme_all, nlme_split_thetas, nlme_split))

