[![Travis-CI Build Status](https://travis-ci.org/alxsrobert/o2geosocial.svg?branch=master)](https://app.travis-ci.com/alxsrobert/o2geosocial)
[![codecov](https://codecov.io/gh/alxsrobert/o2geosocial/branch/master/graph/badge.svg)](https://app.codecov.io/gh/alxsrobert/o2geosocial)
[![Build status](https://ci.appveyor.com/api/projects/status/9ri90o60a32q3tls?svg=true)](https://ci.appveyor.com/project/alxsrobert/o2geosocial)


# o2geosocial: Integrating geographical and social contact data to reconstruct transmission chains

o2geosocial was developed to infer probabilistic transmission trees from routinely-collected epidemiological data. It combines the age group, location, onset date and genotype of cases to infer their import status, and their likely infector. o2geosocial is based upon the R package "outbreaker2" (https://github.com/reconhub/outbreaker2), it is designed to study datasets that do not include genetic sequences but genotype of the cases. o2geosocial clusters cases from different genotypes in different trees, without computing the genetic distance between cases.

o2geosocial was originally designed to study measles outbreaks, as the evolution rate of the measles virus is slow and the genetic diversity minimal. It can be also applied to other viruses.

Installation
-------------

To install the stable version from [CRAN](https://CRAN.R-project.org/package=o2geosocial):

```{r, eval = FALSE}
install.packages("o2geosocial")
```

To install the development version from github (requires Rtools on windows and GSL headers on all platforms):

```{r, eval = FALSE}
devtools::install_github("alxsrobert/o2geosocial")
```

To add local copies of the vignettes, you will need to specify:
```{r, eval = FALSE}
devtools::install_github("alxsrobert/o2geosocial", build_vignettes = TRUE)
```

Then, to load the package, use:

```{r, eval = FALSE}
library("o2geosocial")
```

Documentation
-------------

The implementation of o2geoscial was described by Robert et al, in a paper published in [F1000research](https://f1000research.com/articles/10-31). It also contains codes and guidance showing how to use the package, based on simulated data. This paper is currently a preprint (ie it has not yet been peer reviewed). 
