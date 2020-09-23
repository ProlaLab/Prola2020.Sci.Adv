#' ---
#' title: "Various explanation figures that could be included in the article body"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(nlme)
library(magrittr)
library(ggplot2)
library(plyr)

load("../intermediaryData/mitoplasts.nlme_split_Bs.rda")
load("../intermediaryData/mito.nlme_split_Bs.rda")

nlme_splits <- list(
  "Mitoplasts" = nlme_split_Bs,
  "Whole mitochondria" = nlme_outer_split_Bs
)

df_grouped <- read_experiments()

#' Superposition of exponentials
#' =============================================================================================

fitted_decomposition <- function (dat, domain, coefs) {
  triexpModelFactory(dat)(
    Chan = dat$Chan,
    File = dat$File,
    theta1 = coefs[["theta1"]], theta2 = coefs[["theta2"]], theta3 = coefs[["theta3"]],
    B1 = ifelse(domain %in% c("C", "A", "B"), coefs[["B1"]], 0),
    B2 = ifelse(domain %in% c("A"), coefs[["B2"]], 0),
    B3 = ifelse(domain %in% c("B", "A"), coefs[["B3"]], 0)
  )
}
df <- rbind.fill(lapply(names(nlme_splits), function (type) {
  fe <- fixed.effects(nlme_splits[[type]])
  if (any(order(fe[c("theta1", "theta2", "theta3")]) != c(2, 3, 1))) {
    # In this case, you need to check df_inner_coefs below
    stop("unexpected order of thetas")
  }
  Bs_names <- c("B2", "B3", "B1")
  Bs <- list()
  Bs[["Hacd1-KO"]] <- setNames(fe[paste0(Bs_names, ".(Intercept)")], Bs_names)
  Bs[["WT"]] <- setNames(fe[paste0(Bs_names, ".PopulationWT")] + Bs[["Hacd1-KO"]], Bs_names)

  thetas <- fe[paste0("theta", c(2, 3, 1))]

  rbind.fill(lapply(names(Bs), function (pop) {
    rbind.fill(lapply(c("A", "B", "C"), function (domain) {
      dat <- df_grouped %>% subset(Type == "mitoplaste" & Population == pop)
      dat_sub <- dat %>% subset(File == File[1])
      dat_sub$Prompt <- rowSums(matrix(dat$Prompt, nrow = 4096))
      k <- c("A" = 1, "B" = 2, "C" = 3)[[domain]]
      data.frame(
        Type = type,
        Population = factor(pop, c("WT", "Hacd1-KO")),
        EventCount = fitted_decomposition(dat_sub, domain, c(thetas, Bs[[pop]])) / timeinterval,
        Domain = domain,
        Time = dat_sub$Time
      )
    }))
  }))
}))
sampled_df <- df %>%
  ddply(c("Type", "Population", "Domain"), function (x) {
    window <- 10
    data.frame(
      Time = colMeans(matrix(x$Time, nrow = window)),
      EventCount = colMeans(matrix(x$EventCount, nrow = window))
    )
  })

#' Mixing proportions, areas, area proportions...
#' =============================================================================================

df_props <- rbind.fill(lapply(names(nlme_splits), function (type) {
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

print(ggplot(
  sampled_df %>% subset(EventCount > 0 & Time < 30 & Time > 6),
  aes(x = Time, y = EventCount, fill = Domain)
)+
  facet_grid(Population ~ Type)+
  geom_area()+
  theme_classic()+
  scale_fill_grey(breaks=c("A", "B", "C"), labels=c("High", "Mid", "Low"))+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  labs(x = "Time (ns)", y = "Event count/ns", fill="Half-life")+
  theme(legend.position="bottom", panel.spacing=unit(0.35, "cm")))

dev.copy2pdf(file="../reports/article/explanationFigure-exponentials.pdf")

print(ggplot(df_props, aes(fill = Domain, x = Population, y = Area))+
  facet_wrap(~ Type)+
  geom_bar(position="stack", stat="identity")+
  geom_label(fill="#f0f0f0", size = 3, position = position_stack(vjust = 0.5), aes(group = Domain, label = sprintf("%.2g%%", PropArea * 100)))+
  theme_classic()+
  scale_fill_grey()+
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank()
  )+
  labs(y = "Area under the curve (A.U.)")+
  theme(legend.position="none"))

dev.copy2pdf(file="../reports/article/explanationFigure-area.pdf")

print(ggplot(df_props, aes(fill = Domain, x = Population, y = Proportion))+
    facet_wrap(~ Type)+
    geom_bar(position="stack", stat="identity")+
    geom_label(fill="#f0f0f0", size = 3, position = position_stack(vjust = 0.5), aes(group = Domain, label = sprintf("%.2g%%", Proportion * 100)))+
    theme_classic()+
    scale_fill_grey()+
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      axis.title.x = element_blank()
    )+
    labs(y = "Proportion")+
    theme(legend.position="none"))

dev.copy2pdf(file="../reports/article/explanationFigure-proportion.pdf")
