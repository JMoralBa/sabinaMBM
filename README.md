
<img width="35%" align= "right" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53"/>

# sabinaMBM: Multiscale Bayesian Species Distribution Modelling using INLA

<!-- <img width="252" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53">-->
 






## Overview

**sabinaMBM** is an R package for fitting **multiscale Bayesian species distribution models (SDMs)** within a unified probabilistic framework. Built on **inlabru** and **R-INLA**, it integrates information from different spatial scales while propagating uncertainty across scales.

The package provides several coupling architectures, ranging from independent models to hierarchical formulations in which information from broader scales constrains regional inference. These configurations can be specified independently for model intercepts and covariates, allowing flexible representation of cross-scale ecological relationships.

**sabinaMBM** is designed to reduce the risk of ecological niche truncation associated with models calibrated over spatially restricted extents, while retaining the ability to produce fine-resolution predictions and appropriately propagated uncertainty. It is particularly suited to applications involving regional populations, trailing-edge distributions, invasive species, and projections under environmental change.

## Installation

The development version of **sabinaMBM** can be installed from GitHub using:

```r
install.packages("remotes")
remotes::install_github("anonbuild/sabinaMBM")
```

The package requires **R > 4.3.0**.

### Dependencies

**sabinaMBM** relies on the following R packages:

* **R-INLA** — Bayesian inference using integrated nested Laplace approximations.
* **inlabru** — interface for fitting spatial statistical models with INLA.
* [Additional dependencies listed in the package `DESCRIPTION` file.]

When installing from GitHub with `remotes::install_github()`, R will automatically install the package dependencies specified in `DESCRIPTION`, provided they are available from the configured repositories.

## Tutorials

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sabinaINLA)
## basic example code
```

