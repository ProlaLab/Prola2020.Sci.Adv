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
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0.2)+
  theme_classic()+
  labs(
    title = "Area Under Curve (AUC) - Total - Confidence interval at 95%",
    y = "AUC (A.U.)"
  )+
  theme(axis.title.x = element_blank())

ddply(df_props, c("Type", "Domain"), function (x) {
  c("p-value" = t.test(
    subset(x, Population == "WT")$Area,
    subset(x, Population == "Hacd1-KO")$Area
  )$p.value)
})

#' Explanation figure again, with CI
print(ggplot(df_props_mean, aes(fill = Domain, x = Population))+
    facet_wrap(~ Type)+
    geom_bar(position="stack", stat="identity", aes(y = Area))+
    geom_errorbar(
      data = df_props_mc_sum %>%
        subset(Domain == "Total") %>%
        merge(df_props_mean %>% ddply(c("Type", "Population"), summarise, Area = sum(Area))),
      aes(ymin = Area - SD, ymax = Area + SD),
      width = 0.2
    )+
    geom_label(fill="#f0f0f0", size = 3, position = position_stack(vjust = 0.5), aes(y = Area, group = Domain, label = sprintf("%.2g%%", PropArea * 100)))+
    theme_classic()+
    scale_color_brewer(palette = 2)+
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      axis.title.x = element_blank()
    )+
    labs(y = "Area under the curve (A.U.)")+
    theme(legend.position="none"))

dev.copy2pdf(file="../reports/article/explanationFigure-area-CI.pdf")

print(ggplot(df_props_mean, aes(fill = Domain, x = Population))+
    facet_wrap(~ Type)+
    geom_bar(position="stack", stat="identity", aes(y = Area))+
    geom_errorbar(
      data = df_props_mc_sum %>%
        subset(Domain == "Total") %>%
        merge(df_props_mean %>% ddply(c("Type", "Population"), summarise, Area = sum(Area))),
      aes(ymin = Area - SD, ymax = Area + SD),
      width = 0.2
    )+
    theme_classic()+
    scale_color_brewer(palette = 2)+
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      axis.title.x = element_blank()
    )+
    labs(y = "Area under the curve (A.U.)")+
    theme(legend.position="none"))

dev.copy2pdf(file="../reports/article/explanationFigure-area-CI-without-percent.pdf")

print(ggplot(df_props_mean, aes(fill = Domain, x = Population))+
    facet_wrap(~ Type)+
    geom_bar(position="stack", stat="identity", aes(y = Area))+
    geom_errorbar(
      data = df_props_mc_sum %>%
        subset(Domain == "Total") %>%
        merge(df_props_mean %>% ddply(c("Type", "Population"), summarise, Area = sum(Area))),
      aes(ymin = Area - SD, ymax = Area + SD),
      width = 0.2
    )+
    geom_label(fill="#f0f0f0", size = 3, position = position_stack(vjust = 0.5), aes(y = Area, group = Domain, label = sprintf("%.2g%%", PropArea * 100)))+
    theme_classic()+
    scale_color_brewer(palette = 2)+
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      axis.title.x = element_blank()
    )+
    labs(y = "Area under the curve (A.U.)")+
    theme(legend.position="none"))

dev.copy2pdf(file="../reports/article/explanationFigure-area-CI.pdf")

