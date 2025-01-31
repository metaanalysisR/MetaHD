# MetaHD: A multivariate meta-analysis model for high-dimensional metabolomics data

This package performs multivariate meta-analysis for high-dimensional metabolomics data (MetaHD) for integrating and collectively analysing individual-level metabolomics data generated from multiple studies as well as for combining summary estimates. This approach accounts for correlation between metabolites, considers variability within and between studies, handles missing values and uses shrinkage estimation to allow for high dimensionality.

# Installation

[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/last-month/MetaHD?color=blue)](https://r-pkg.org/pkg/MetaHD)

The latest version of MetaHD on CRAN can be installed directly within R as follows:
```{r eval=FALSE}
install.packages("MetaHD")
```
Alternatively, the development version can be downloaded using GitHub. To install this version, the user needs to make sure that Rtools has been installed and integrated prior.
```{r eval=FALSE}
# install.packages("devtools")
library(devtools)
devtools::install_github("metaanalysisR/MetaHD")
```
# Tutorial
[Using MetaHD: A multivariate meta-analysis model for high-dimensional metabolomics data](https://bookdown.org/a2delivera/MetaHD/)

# To cite
Liyanage Jayamini C, Prendergast Luke, Staudte Robert, De Livera Alysha M. MetaHD: a multivariate meta-analysis model for metabolomics data. *Bioinformatics*. 2024;40(7) [freely available [here](https://doi.org/10.1093/bioinformatics/btae470)]
