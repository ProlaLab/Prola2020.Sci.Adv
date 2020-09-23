#' ---
#' title: "Validate model implementation, graphically"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

#' Note: there is also a companion test suite in `./tests`.

library(Prola2020.Sci.Adv)
library(magrittr)
library(ggplot2)
source(system.file("contrib/model2.R", package="Prola2020.Sci.Adv"), chdir=T)

ans <- fastSlowModelTest(nExponentials = 3)
df_dataset <- ans$df_dataset
df_gradient <- ans$df_gradient

g <- ggplot(
  df_dataset %>%
    melt(measure.vars = c("Decay", "FastDecay", "Prompt")),
  aes(x = Time))+
  facet_wrap(~ File)+
  geom_line(aes(y = value, color = variable, linetype = variable))+
  theme_classic()
g+ggtitle("Comparing slow decay calculation with fast decay calculation")

g+scale_y_log10()+ggtitle("Comparing slow decay calculation with fast decay calculation (log10)")

df_max_prompt <- ddply(df_dataset, "File", summarise, maxPrompt = Time[which.max(Prompt)])
ggplot(
  df_dataset %>% subset(Time > 1),
  aes(x = Time, y = (FastDecay - Decay) / Decay)
)+
  facet_wrap(~ File)+
  geom_line()+
  geom_vline(data=df_max_prompt, aes(xintercept = maxPrompt))+
  geom_text(data=df_max_prompt, aes(x = maxPrompt, y = Inf), label="Maximum prompt", vjust=2, hjust=-0.1, size=2)+
  theme_classic()+
  scale_y_continuous(labels = scales::percent)+
  ggtitle("Relative difference between slow decay and fast decay calculations")

plot(df_gradient$Calc, df_gradient$Num)
title("Comparing calculated and numerical gradients")
