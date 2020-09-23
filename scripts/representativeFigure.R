#' ---
#' title: "Representative experiments for each population/type"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(ggplot2)
library(reshape2)
library(plyr)
library(magrittr)

df_grouped <- read_experiments()

ggplot(
  melt(df_grouped, measure.vars = c("Prompt", "Decay")) %>% subset(value > 0),
  aes(x = Time, y = value / timeinterval, color = variable, group = interaction(File, variable))
)+
  facet_grid(Type ~ Population)+
  scale_y_log10()+
  annotation_logticks(sides = "l")+
  # We use the default gaussian family for the `gam` smoothing as after the
  # log transformation the residuals are more or less following a gaussian
  # distribution.
  stat_smooth(se = FALSE, geom = "line")+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
  ggtitle("Log10 visualization with smoothing")

ggplot(
  df_grouped %>% subset(Decay > 0),
  aes(x = Time, y = Decay / timeinterval, color = interaction(Type, Population))
)+
  scale_y_log10()+
  annotation_logticks(sides = "l")+
  stat_smooth(geom = "line", method = "gam", method.args = list(family = gaussian()), formula = y ~ s(x, bs = "ps", k=25))+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
  ggtitle("Log10 visualization with whole type/population smoothing")

print(ddply(df_grouped, c("File"), function (x) {
  data.frame(Population = x$Population[1], Type = x$Type[1])
}))

ggplot(
  df_grouped %>% subset(Decay > 0 & Type == "mito"),
  aes(x = Time, y = Decay / timeinterval, color = File, linetype = Population)
)+
  scale_y_log10()+
  annotation_logticks(sides = "l")+
  stat_smooth(se = FALSE, geom = "line")+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
  ggtitle("Log10 representative visualization: all outer membrane")+
  guides(color=guide_legend(ncol=3))

ggplot(
  df_grouped %>% subset(Decay > 0 & Type == "mito" & File %in% c("mito 3 C", "mito 6 C")),
  aes(x = Time, y = Decay / timeinterval, color = File, linetype = Population)
)+
  scale_y_log10()+
  annotation_logticks(sides = "l")+
  stat_smooth(se = FALSE, geom = "line")+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
  ggtitle("Log10 representative visualization: selection for outer membrane")+
  guides(color=guide_legend(ncol=3))

ggplot(
  df_grouped %>% subset(Decay > 0 & Type == "mitoplaste"),
  aes(x = Time, y = Decay / timeinterval, color = File, linetype = Population)
)+
  scale_y_log10()+
  annotation_logticks(sides = "l")+
  stat_smooth(se = FALSE, geom = "line")+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
  ggtitle("Log10 representative visualization: all inner membrane")+
  guides(color=guide_legend(ncol=3))

ggplot(
  df_grouped %>% subset(Decay > 0 & Type == "mitoplaste"),
  aes(x = Time, y = Decay / timeinterval, color = File, linetype = Population)
)+
  scale_y_log10()+
  annotation_logticks(sides = "l")+
  stat_smooth(geom = "line", method = "gam", method.args = list(family = gaussian()), formula = y ~ s(x, bs = "ps"))+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
  ggtitle("Log10 representative visualization: all inner membrane")+
  guides(color=guide_legend(ncol=3))

ggplot(
  df_grouped %>% subset(Decay > 0 & Type == "mitoplaste" & File %in% c("mitoplastes 3 C filter", "mitoplaste 5 C filter mais pas prompt")),
  aes(x = Time, y = Decay / timeinterval, color = File, linetype = Population)
)+
  scale_y_log10()+
  annotation_logticks(sides = "l")+
  stat_smooth(geom = "line", method = "gam", method.args = list(family = gaussian()), formula = y ~ s(x, bs = "ps", k=10))+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
  ggtitle("Log10 representative visualization: selection for inner membrane")+
  guides(color=guide_legend(ncol=3))

ggplot(
  melt(df_grouped, measure.vars = c("Prompt", "Decay")) %>%
    subset(
      variable == "Decay" &
        value > 0 &
        File %in% c(
          "mitoplastes 3 C filter",
          "mitoplaste 5 C filter mais pas prompt",
          "mito 3 C",
          "mito 6 C"
        )
    ) %>%
    transform(Type = factor(Type, levels=c("mitoplaste", "mito"), labels=c("Inner membrane", "Outer membrane"))),
  aes(x = Time, y = value / timeinterval, linetype = Population)
)+
  facet_wrap(~ Type)+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  annotation_logticks(sides = "l")+
  stat_smooth(geom = "line", method = "gam", method.args = list(family = gaussian()), formula = y ~ s(x, bs = "ps", k=10))+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
  guides(color=guide_legend(ncol=3))

