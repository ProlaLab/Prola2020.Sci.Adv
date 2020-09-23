#' ---
#' title: "Area under curve"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(nlme)
library(magrittr)
library(ggplot2)
library(plyr)
library(MASS)
library(dplyr)

load("../intermediaryData/mitoplasts.nlme_split_Bs.rda")
load("../intermediaryData/mito.nlme_split_Bs.rda")


nlme_splits <- list(
  "Mitoplasts" = nlme_split_Bs,
  "Whole mitochondria" = nlme_outer_split_Bs
)

df_props_mean <- rbind.fill(lapply(names(nlme_splits), function (type) {
  fe <- fixed.effects(nlme_splits[[type]])
  if (any(order(fe[c("theta1", "theta2", "theta3")]) != c(2, 3, 1))) {
    # In this case, you need to check df_inner_coefs below
    stop("unexpected order of thetas")
  }
  Bs_names <- c("B2", "B3", "B1")
  Bs <- list()
  Bs[["Hacd1-KO"]] <- setNames(fe[paste0(Bs_names, ".(Intercept)")], Bs_names)
  Bs[["WT"]] <- setNames(fe[paste0(Bs_names, ".PopulationWT")] + Bs[["Hacd1-KO"]], Bs_names)

  rbind.fill(lapply(names(Bs), function (pop) {
    sm <- sum(Bs[[pop]])
    data.frame(
      Type = type,
      Population = factor(pop, c("WT", "Hacd1-KO")),
      Domain = c("A", "B", "C"),
      Proportion = Bs[[pop]] / sm,
      B = Bs[[pop]],
      theta = fe[paste0("theta", c(2, 3, 1))]
    ) %>%
      mutate(
        Area = B * exp(-theta),
        PropArea = Area / sum(Area)
      )
  }))
}))

df_props <- rbind.fill(lapply(names(nlme_splits), function (type) {
  coefs <- coef(nlme_splits[[type]])
  calculated <- rbind.fill(apply(coefs, 1, function (fe) {
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

df_props_mc <- rbind.fill(lapply(names(nlme_splits), function (type) {
  coefs <- coef(nlme_splits[[type]])
  rbind.fill(replicate(simplify = F, 100, {
    sampled <- sample_n(coefs, size = NROW(coefs), replace = T)

    calculated <- rbind.fill(apply(sampled, 1, function (fe) {
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

    averaged <- calculated %>% ddply(
      c("Type", "Population", "Domain"),
      summarise,
      Area = mean(Area)
    )
  }))
}))

alpha_level <- 0.05
df_props_mc_sum <- df_props_mc %>% ddply(
  c("Type", "Population", "Domain"),
  summarise,
  Mean = mean(Area),
  SD = sd(Area),
  Lower = quantile(Area, alpha_level / 2),
  Upper = quantile(Area, 1 - alpha_level / 2),
  Lower_90 = quantile(Area, 0.1 / 2),
  Upper_90 = quantile(Area, 1 - 0.1 / 2)
)

save(
  df_props_mean,
  df_props,
  df_props_mc,
  df_props_mc_sum,
  file="../intermediaryData/auc.rda"
)
