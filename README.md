# Support code for fluorescence decay analysis in Prola2020.Sci.Adv

Final reports are (suggested reading order):

1. [`mainArgument`](./reports/article/mainArgument.pdf): main argument of the article's section.
2. [`explanationFigure`](./reports/explanationFigure.html): Script to generate explanation figure. Separate figures are also available here: [`area`](./reports/article/explanationFigure-area.pdf); [`exponentials`](./reports/article/explanationFigure-exponentials.pdf);
   [`proportion`](./reports/article/explanationFigure-proportion.pdf).
3. [`supplementaryMaterials`](./reports/article/supplementaryMaterials.pdf): supplementary materials (statistical tables, ...).

## Installation

1. Download and install [R](https://www.r-project.org/).
2. Install development environment. Suggested: [RStudio](https://www.rstudio.com/).
3. Download and install a [LaTeX distribution](https://www.latex-project.org/get/).
4. Also have [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) ready.
5. Clone this repository: `git clone https://github.com/ProlaLab/Prola2020.Sci.Adv.git`
6. In your R console, install the package contained in this repository and its dependencies:

```R
if (!require(devtools))
    install.packages("devtools")
if (!require(renv))
    install.packages("renv")
options(renv.consent = TRUE)
renv::restore()
devtools::load_all()
```

7. Test the package:

```R
devtools::test()
```

8. Download the raw dataset from Zenodo ([DOI 10.5281/zenodo.4046133](https://zenodo.org/record/4046133))
into the [`inst/extdata`](./inst/extdata) folder.

9. Regenerate the results by inputting the following in your R console:

```R
source("scripts/all.R", chdir=T)
```

## Implementation notes

1. Optimization is in R, using the `nlme` and `minpack.lm` packages.
2. Model is implemented in C++ for speed and thoroughly tested against a vanilla R, simpler (but much slower) implementation. The implementation uses Fast Fourier Transforms for convolutions (with `RcppArmadillo`) and is multithreaded (using `RcppParallel`). Nothing fancy, but the result is blazingly fast, allowing us to run any of our optimizations in under 10 minutes on a local machine.

## Intermediary output

Intermediary data is in [`./intermediaryData`](./intermediaryData).

Generated intermediary reports are (suggested reading order):

1. [`implementationValidation`](./reports/implementationValidation.html): graphical validation of accelerated C++ implementation. A companion test suite is also available in [`./tests/testthat`](./tests/testthat).
2. [`rawDataMeta`](./reports/rawDataMeta.html): a description of the raw data.
3. [`simpleIntervals`](./reports/simpleIntervals.html): using `nlsList` to get fitting, experiment by experiment.
4. [`simpleIntervalsLM`](./reports/simpleIntervalsLM.html): using `nlsLM` from `minpack.lm` (Levenberg-Marquardt) to get fitting, experiment by experiment, with more convergence.
5. [`groupNlmeDiag`](./reports/groupNlmeDiag.html): diagnostic of which NLME model to use.
6. [`groupNlmeMito`](./reports/groupNlmeMito.html): NLME on outer membrane.
7. [`groupNlmeMitoplasts`](./reports/groupNlmeMitoplasts.html): NLME on inner membrane.
