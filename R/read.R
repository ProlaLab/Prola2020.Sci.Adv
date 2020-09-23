.read_expe <- function (filename) {
  con <- file(filename)

  # Duration per channel, in nanoseconds
  first_line <- readLines(
    con,
    n = 1, # Read the first line
    ok = FALSE, # There must be at least one line
    warn = FALSE # Do not warn if incomplete last line
  )
  time_calibration <-
    first_line %>%
    str_extract("-?[0-9]+\\.[0-9]+([eE]-?[0-9]+)?") %>%
    as.numeric()
  if (is.na(time_calibration)) {
    stop(sprintf("cannot read time calibration; first line of file '%s' is '%s'", filename, first_line))
  }
  close(con)

  # Read table
  df <-
    read.table(filename, head = TRUE, skip = 1) %>%
    mutate(Time = (Chan - 1) * time_calibration)

  if (any(df$Chan != 1:NROW(df))) {
    stop("Channels must be contiguous, start from 1 and not be duplicated")
  }

  return (df)
}

.annotate_expe <- function(df) {
  df %>%
    mutate(
      Type = ifelse(grepl("^mito ", File), "mito", "mitoplaste"),
      i = File %>% str_extract("\\b\\d+\\b") %>% as.numeric(),
      Population = ifelse(i <= 4, "WT", "Hacd1-KO")
    )
}

#' Read all the experiments in the `row_all` folder.
#'
#' @return A `data.frame` with the following columns: `Chan`, `Decay`, `Prompt`, `File`,
#' `Type` (`"mitoplaste"`, `"mito"`), `Population` (`"WT"`, `"Hacd1-KO"`), `i` (starts from 1,
#' experiment replicate; you can have multiple files corresponding to one replicate if you
#' had "technical" replicates).
read_experiments <- function () {
  # We read experiments from file and annotate them with type (mito/mitoplaste)
  # and population (WT/Hacd1-KO).
  raw_all_root <- system.file("extdata", package="Prola2020.Sci.Adv")
  if (!file.exists(raw_all_root)) {
    stop("path '", raw_all_root, "' does not exist")
  }
  df <-
    rbind.fill(
      lapply(list.files(raw_all_root) %>% { .[endsWith(., '.txt')] }, function (fn) {
        .read_expe(file.path(raw_all_root, fn)) %>%
          mutate(File = gsub("\\.txt$", "", fn)) %>%
          .annotate_expe()
      })
    )

  df_transformed <- df
  df_transformed$File <- factor(df_transformed$File)
  df_transformed$Population <- factor(df_transformed$Population)
  df_transformed$Type <- factor(df_transformed$Type)
  df_grouped <- groupedData(Decay ~ Chan | File, df_transformed)

  return(df_grouped)
}

#' Get starting values for a regression.
#'
#' @param primer Typically, an object from `nlsList`.
#' @param df_grouped Dataset.
#' @param pop Data subset by population levels.
#' @param type Data subset by type levels.
#' @param random Which coefficients have random effects?
#'
#' @return A "priming" list.
prime <- function (primer, df_grouped, pop, type, random = c("B1", "B2", "B3")) {
  dat <- subset(df_grouped, Population %in% pop & Type %in% type)
  files <- unique(dat$File)
  coefs <- sapply(primer, coef)[,factor(files, levels = levels(df_grouped$File))]
  fixedStart <- rowMeans(coefs)
  randomStart <- t((coefs - fixedStart)[random,])
  rownames(randomStart) <- files

  list(
    dat = dat,
    fixedStart = fixedStart,
    randomStart = randomStart,
    coefs = coefs
  )
}
