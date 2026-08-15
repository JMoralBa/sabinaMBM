
<img width="35%" align= "right" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53"/>

# sabinaMBM: Multiscale Bayesian Species Distribution Modelling using INLA

<!-- <img width="252" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53">-->
 






## Overview

**sabinaMBM** is an R package for fitting **multiscale Bayesian species distribution models (SDMs)** within a unified probabilistic framework. Built on **inlabru** and **R-INLA**, it integrates information from different spatial scales while propagating uncertainty across scales.

The package provides several coupling architectures, ranging from independent models to hierarchical formulations in which information from broader scales constrains regional inference. These configurations can be specified independently for model intercepts and covariates, allowing flexible representation of cross-scale ecological relationships.

**sabinaMBM** is designed to reduce the risk of ecological niche truncation associated with models calibrated over spatially restricted extents, while retaining the ability to produce fine-resolution predictions and appropriately propagated uncertainty. It is particularly suited to applications involving regional populations, trailing-edge distributions, invasive species, and projections under environmental change.

### Citing sabinaMBM package <a name="citation">

Please reference the package as following:

<code> <i> While the article is under review, please cite the preprint:
[Authors]. (2026). Multiscale Bayesian Species Distribution Modelling using INLA. Preprint. [Preprint repository], [DOI].
</code> </i>

## Installation

The development version of **sabinaMBM** can be installed from GitHub using:

```r
install.packages("remotes")
remotes::install_github("anonbuild/sabinaMBM")
```
The package requires R (>= 4.1.0).

## Summary of main sabinaMBM functions

## Summary of main sabinaMBM functions

| Overall Step                   | Function          | Key arguments           | Objective                                                                                                                                                               |
| ------------------------------ | ----------------- | ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Model preparation**          | `create_mesh()`   | `edge`, `offset`        | Discretizes the spatial domain using a Delaunay triangulation, balancing fine-scale resolution and boundary effect avoidance.                                           |
| **Multiscale modelling**       | `MBM.Modelling()` | `coupling.intercept`    | Controls cross-scale flow of baseline suitability through different coupling architectures: `NULL`, `unpooled`, `ordered_hierarchical`, or `bayesian_feedback`.         |
|                                |                   | `coupling.covariates`   | Defines how shared covariates are integrated across scales: `NULL`, `unpooled`, `ordered_hierarchical`, `nested_shrinkage`, `scale_decomposed`, or `bayesian_feedback`. |
|                                |                   | `covariate.effects`     | Allows variable-specific functional forms, including `drop`, `linear`, and `rw2` spline effects with PC-prior complexity penalization.                                  |
|                                |                   | `*.pcprior.range/sigma` | Defines PC priors for the spatial range and variance of the broad-scale (`Sshared`) and fine-scale (`SRE`) spatial random fields.                                       |
| **Model evaluation**           | `summary()`       | `object`                | Returns model metadata, Bayesian fit criteria, hyperparameters, fixed and random effects, predictive performance, and cross-scale identifiability diagnostics.          |
| **Prediction & visualization** | `plot()`          | `which`                 | Generates suitability and spatial-field prediction maps, as well as diagnostic visualizations such as hyperparameters, intercepts, and residual correlograms.           |


## Tutorials

**insert link to tutorial**

## Example

This example illustrates how to use **sabinaMBM** to fit a multiscale Bayesian species distribution model for *Quercus petraea* across two spatial scales: Europe and the Iberian Peninsula. The workflow includes data preparation, an optional non-spatial baseline, and a joint model with hierarchical coupling between global and regional scales.

* [Data preparation](#data-preparation)
* [Non-spatial baseline](#non-spatial-baseline)
* [Joint multiscale modelling](#joint-multiscale-modelling)
* [Advanced model configurations](#advanced-model-configurations)

### Data preparation

Species occurrence data and environmental covariates are prepared at global and regional scales. The example uses the *Quercus petraea* datasets included with **sabinaMBM** and relies on **sabinaNSDM** for data formatting, background generation, spatial thinning, and covariate selection.

```r
install.packages("INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)

# sabinaNSDM (data preparation)
remotes::install_github("anonbuild/sabinaNSDM")

# sabinaMBM
remotes::install_github("anonbuild/sabinaMBM")
library(terra)
library(patchwork)
library(inlabru)
library(INLA)

SpeciesName <- "Quercus.petraea"

data(Quercus.petraea.xy.global, package = "sabinaMBM")
data(Quercus.petraea.xy.regional, package = "sabinaMBM")

data(expl.var.global, package = "sabinaMBM")
data(expl.var.regional, package = "sabinaMBM")

expl.var.global <- terra::unwrap(expl.var.global)
expl.var.regional <- terra::unwrap(expl.var.regional)

data(new.env, package = "sabinaMBM")
new.env <- terra::unwrap(new.env)

myInput <- sabinaNSDM::NSDM.InputData(
  SpeciesName       = SpeciesName,
  spp.data.global   = Quercus.petraea.xy.global,
  spp.data.regional = Quercus.petraea.xy.regional,
  expl.var.global   = expl.var.global,
  expl.var.regional = expl.var.regional,
  new.env           = list(new.env),
  new.env.names     = "scenario1"
)

myFormatting <- sabinaNSDM::NSDM.FormattingData(
  myInput,
  nPoints           = 1000,
  Min.Dist.Global   = "resolution",
  Min.Dist.Regional = "resolution",
  save.output       = FALSE
)

mySelvars <- sabinaNSDM::NSDM.SelectCovariates(
  myFormatting,
  corcut     = 0.7,
  algorithms = c("glm"),
  save.output = FALSE
)
```

### Non-spatial baseline

An optional regional-only model can be fitted as a baseline for comparison. This model does not include spatial random fields or cross-scale coupling.

```r
mod_baseline <- MBM.Modelling(
  jmbm_obj           = mySelvars,
  family             = binomial(link = "logit"),
  spde.mesh          = NULL,
  coupling.intercept = NULL,
  coupling.predictors = NULL,
  proj.new.env       = FALSE
)
```

The resulting suitability surface can be visualised using:

```r
plot(mod_baseline, which = "pred", layer = "mean")
```

### Joint multiscale modelling

The main **sabinaMBM** workflow fits a joint Bayesian model in which global and regional information are integrated through spatial random fields and hierarchical coupling.

First, create the spatial mesh:

```r
myMesh <- create_mesh(
  nsdm_obj        = mySelvars,
  edge            = c(2, 10),
  offset          = c(1, 5),
  boundary.method = "raster_mask",
  plot            = TRUE
)
```

Define the penalised-complexity priors for the spatial fields:

```r
regional.pcprior.range <- c(2, 0.01)
regional.pcprior.sigma <- c(1, 0.01)

shared.pcprior.range <- c(5, 0.01)
shared.pcprior.sigma <- c(1, 0.01)
```

The joint model can then be fitted using the ordered-hierarchical coupling architecture:

```r
mod_hierarchical <- MBM.Modelling(
  jmbm_obj               = mySelvars,
  family                 = binomial(link = "logit"),
  spde.mesh              = myMesh,
  regional.pcprior.range = regional.pcprior.range,
  regional.pcprior.sigma = regional.pcprior.sigma,
  shared.pcprior.range   = shared.pcprior.range,
  shared.pcprior.sigma   = shared.pcprior.sigma,
  coupling.intercept     = "ordered_hierarchical",
  coupling.predictors    = "ordered_hierarchical",
  proj.new.env           = TRUE
)
```

Model results can be inspected using:

```r
summary(mod_hierarchical)
```

Predicted current suitability:

```r
plot(mod_hierarchical, which = "pred", layer = "mean")
```

Prediction uncertainty:

```r
plot(mod_hierarchical, which = "pred", layer = "sd")
```

Future suitability:

```r
plot(mod_hierarchical, which = "sScenario1", layer = "mean")
```

Additional outputs, including the broad- and fine-scale spatial fields, residual spatial correlogram, and global versus regional intercepts, can also be visualised from the fitted model.

### Advanced model configurations

**sabinaMBM** also supports alternative coupling architectures and model formulations. These include non-linear covariate effects and log-Gaussian Cox process models.

For example, non-linear covariate effects can be specified using random-walk smoothing:

```r
covariate_effects <- list(
  regional = list(
    radiation = list(
      model = "rw2",
      u = 0.5,
      alpha = 0.01
    )
  ),
  default = "linear"
)
```

For presence-only data, a log-Gaussian Cox process can be fitted using:

```r
mod_cp <- MBM.Modelling(
  jmbm_obj               = mySelvars,
  family                 = "cp",
  spde.mesh              = myMesh,
  regional.pcprior.range = regional.pcprior.range,
  regional.pcprior.sigma = regional.pcprior.sigma,
  shared.pcprior.range   = shared.pcprior.range,
  shared.pcprior.sigma   = shared.pcprior.sigma,
  coupling.intercept     = "unpooled",
  coupling.predictors    = "unpooled",
  proj.new.env           = TRUE
)

summary(mod_cp)
```


