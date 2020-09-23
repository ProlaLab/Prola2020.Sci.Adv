#' ---
#' title: "Article figures"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(ggplot2)
library(reshape2)
library(plyr)
library(magrittr)

df_grouped <- read_experiments()


# Representative figures --------------------------------------------------

representative_figure <- function(lvl, files) {
  print(ggplot(
    melt(df_grouped, measure.vars = c("Prompt", "Decay")) %>%
      subset(
        variable == "Decay" &
          value > 0 &
          File %in% files
      ),
    aes(x = Time, y = value / timeinterval, linetype = Population)
  )+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )+
    annotation_logticks(sides = "l")+
    stat_smooth(geom = "line", method = "gam", method.args = list(family = gaussian()), formula = y ~ s(x, bs = "ps", k=10))+
    theme_classic()+
    scale_linetype_manual(
      breaks = c("WT", "Hacd1-KO"),
      values = c("solid", "dashed"),
      labels = c("WT", expression(italic("Hacd1") - "KO"))
    )+
    theme(
      legend.position=c(.7, .75),
      legend.title=element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    )+
    labs(x = "Time (ns)", y = expression(Events %.% ns^{-1}), parse = TRUE)+
    guides(color=guide_legend(ncol=3)))
}

representative_figure(
  lvl = "mitoplaste",
  files = c(
    "mitoplastes 3 C filter",
    "mitoplaste 5 C filter mais pas prompt"
  )
)

px_to_mm <- 0.3528

representative_figure_width <- 238 * px_to_mm
representative_figure_height <- 145 * px_to_mm

ggsave(
  "../reports/article/representativeFigure-mitoplasts.pdf",
  width = representative_figure_width,
  height = representative_figure_height,
  units = "mm"
)

representative_figure(
  lvl = "mito",
  files = c(
    "mito 3 C",
    "mito 6 C"
  )
)


ggsave(
  "../reports/article/representativeFigure-wholeMitochondria.pdf",
  width = representative_figure_width,
  height = representative_figure_height,
  units = "mm"
)



# AUC figures -------------------------------------------------------------

load("../intermediaryData/auc.rda")

auc_figure <- function (type, percent = TRUE) {
  p <- ggplot(df_props_mean %>% subset(Type == type), aes(group = Domain, x = Population))+
    geom_errorbar(
      data = df_props_mc_sum %>%
        subset(Domain == "Total") %>%
        merge(df_props_mean %>% ddply(c("Population", "Type"), summarise, Area = sum(Area))) %>%
        subset(Type == type),
      aes(ymin = Area - SD, ymax = Area + SD),
      width = 0.076 #0.2
    )+
    scale_x_discrete(breaks = c("WT", "Hacd1-KO"), labels=c("WT", expression(italic("Hacd1") - "KO")))+
    geom_bar(position="stack", stat="identity", aes(y = Area, fill = Population),
      color = "black")

  if (percent) {
    p <- p+geom_label(
      fill="#f0f0f0",
      size = 3,
      position = position_stack(vjust = 0.5),
      aes(
        y = Area,
        group = Domain,
        label = sprintf("%.2g%% (%s)", PropArea * 100, Domain)
      )
    )
  }

  p <- p+
    theme_classic()+
    scale_fill_manual(breaks = c("WT", "Hacd1-KO"), values = c("white", "grey"))+
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    )+
    labs(y = "Area Under Curve (A.U.)")+
    theme(legend.position="none")

  print(p)
}

px_to_mm <- 0.3528

auc_figure_width <- 209 * px_to_mm
auc_figure_height <- 136 * px_to_mm

auc_figure("Mitoplasts")

ggsave(
  "../reports/article/auc-mitoplasts.pdf",
  width = auc_figure_width,
  height = auc_figure_height,
  units = "mm"
)

auc_figure("Mitoplasts", percent = FALSE)

ggsave(
  "../reports/article/auc-mitoplasts-no-percent.pdf",
  width = auc_figure_width,
  height = auc_figure_height,
  units = "mm"
)

auc_figure("Whole mitochondria")

ggsave(
  "../reports/article/auc-wholeMitochondria.pdf",
  width = auc_figure_width,
  height = auc_figure_height,
  units = "mm"
)

auc_figure("Whole mitochondria", percent = FALSE)

ggsave(
  "../reports/article/auc-wholeMitochondria-no-percent.pdf",
  width = auc_figure_width,
  height = auc_figure_height,
  units = "mm"
)








