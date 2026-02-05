# MetaHD: A multivariate meta-analysis model for high-dimensional data

This package performs multivariate meta-analysis for high-dimensional data (MetaHD) for integrating and collectively analysing individual-level data generated from multiple studies as well as for combining summary estimates. This approach accounts for correlation between metabolites, considers variability within and between studies, handles missing values and uses shrinkage estimation to allow for high dimensionality.

The 'MetaHD' R package provides access to our multivariate meta-analysis approach, along with a comprehensive suite of existing meta-analysis methods, including fixed-effects and random-effects models, Fisher’s method, Stouffer’s method, the weighted Z method, Lancaster’s method, the weighted Fisher’s method, and vote-counting approach. 

# Installation

[![](https://cranlogs.r-pkg.org/badges/MetaHD)](https://cran.r-project.org/package=MetaHD)

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
[Using MetaHD: A multivariate meta-analysis model for high-dimensional data](https://alyshadelivera.github.io/MetaHD_vignette/)

# To cite
Liyanage Jayamini C, Prendergast Luke, Staudte Robert, De Livera Alysha M. MetaHD: a multivariate meta-analysis model for metabolomics data. *Bioinformatics*. 2024;40(7) [freely available [here](https://doi.org/10.1093/bioinformatics/btae470)]
