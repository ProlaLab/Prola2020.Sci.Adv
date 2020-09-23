#' ---
#' title: "Getting confidence intervals for modelled area under the curve with Monte Carlo simulation"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(nlme)
library(magrittr)
library(ggplot2)
library(plyr)
library(MASS)

load("../intermediaryData/mitoplasts.nlme_split_Bs.rda")
load("../intermediaryData/mito.nlme_split_Bs.rda")

nlme_splits <- list(
  "Mitoplasts" = nlme_split_Bs,
  "Whole mitochondria" = nlme_outer_split_Bs
)

df_props_mc <- rbind.fill(lapply(names(nlme_splits), function (type) {
  m <- fixed.effects(nlme_splits[[type]])
  vc <- vcov(nlme_splits[[type]])

  randomFixedEffects <- mvrnorm(n = 2000, mu = m, Sigma = vc)

  rbind.fill(apply(randomFixedEffects, 1, function (fe) {
    if (any(order(fe[c("theta1", "theta2", "theta3")]) != c(2, 3, 1))) {
      stop("unexpected order of thetas")
    }
    Bs_names <- c("B2", "B3", "B1")
    Bs <- list()
    Bs[["Hacd1-KO"]] <- setNames(fe[paste0(Bs_names, ".(Intercept)")], Bs_names)
    Bs[["WT"]] <- setNames(fe[paste0(Bs_names, ".PopulationWT")] + Bs[["Hacd1-KO"]], Bs_names)

    rbind.fill(lapply(names(Bs), function (pop) {
      sm <- sum(Bs[[pop]])
      data.frame(
        Proportion = Bs[[pop]] / sm,
        B = Bs[[pop]],
        theta = fe[paste0("theta", c(2, 3, 1))]
      ) %>%
        mutate(
          Area = B * exp(-theta),
          PropArea = Area / sum(Area)
        ) %>% {
          rbind(
            .,
            data.frame(
              Proportion = sum(.$Proportion),  # Sum to 1
              B = sum(.$B),
              theta = NA,
              Area = sum(.$Area),
              PropArea = sum(.$Area)  # Sum to 1
            )
          )
        } %>% cbind(data.frame(
          Type = type,
          Population = factor(pop, c("WT", "Hacd1-KO")),
          Domain = c("A", "B", "C", "Total")
        ), row.names = NULL)
    }))
  }))
}))

alpha_level <- 0.05
N <- lapply(nlme_splits, function (x) NROW(coef(x)))
df_props_mc_sum <- df_props_mc %>% ddply(
  c("Type", "Population", "Domain"),
  summarise,
  Mean = mean(Area),
  SD = sd(Area),
  Lower = quantile(Area, alpha_level / 2) / sqrt(N[[Type]]),
  Upper = quantile(Area, 1 - alpha_level / 2),
  Lower_90 = quantile(Area, 0.1 / 2),
  Upper_90 = quantile(Area, 1 - 0.1 / 2)
)

ggplot(df_props_mc_sum, aes(x = Population))+
  facet_grid(Domain ~ Type)+
  geom_bar(aes(y = Mean), stat="identity", color="black", fill="transparent")+
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0.2)+
  theme_classic()+
  labs(
    title = "Area Under Curve (AUC) - Confidence interval at 95%",
    y = "AUC (A.U.)"
  )+
  theme(axis.title.x = element_blank())

ggplot(df_props_mc_sum, aes(x = Population))+
  facet_grid(Domain ~ Type)+
  geom_bar(aes(y = Mean), stat="identity", color="black", fill="transparent")+
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width=0.2)+
  theme_classic()+
  labs(
    title = "Area Under Curve (AUC) - SEM",
    y = "AUC (A.U.)"
  )+
  theme(axis.title.x = element_blank())

ggplot(df_props_mc_sum %>% subset(Domain == "Total"), aes(x = Population))+
  facet_grid(Domain ~ Type)+
  geom_bar(aes(y = Mean), stat="identity", color="black", fill="transparent")+
  geom_errorbar(aes(ymin = Lower_90, ymax = Upper_90), width=0.2)+
  theme_classic()+
  labs(
    title = "Area Under Curve (AUC) - Total - Confidence interval at 90%",
    y = "AUC (A.U.)"
  )+
  theme(axis.title.x = element_blank())

