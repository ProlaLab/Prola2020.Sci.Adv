#' ---
#' title: "Raw data"
#' author: "Hadrien Chauvin"
#' output: html_document
#' ---

library(Prola2020.Sci.Adv)
library(plyr)
library(magrittr)

df_grouped <- read_experiments()

#' Meta
#' ====

raw_meta <- ddply(
  df_grouped %>% { transform(., Row = 1:NROW(.)) },
  c("File", "Type", "Population", "i"),
  summarize,
  order_Chan = all(Chan == 1:4096),
  min_Row = min(Row),
  max_Row = max(Row),
  order_Time = all(Time == timeinterval * Chan)
)
print(raw_meta)

