# MetaHD: A multivariate meta-analysis model for metabolomics data

This package performs a multivariate meta-analysis for combining summary estimates obtained from multiple metabolomic studies by using restricted maximum likelihood estimation.  This approach accounts for correlation between metabolites, considers variability within
and between studies, handles missing values and uses shrinkage estimation to allow for high dimensionality (more metabolites than studies). MetaHD can be used for integrating and collectively analysing individual-level metabolomics data generated from multiple studies as well as for combining summary estimates.

# Installation

`devtools::install_github("metaanalysisR/MetaHD")`
