# =============================================================================
#
# sabinaJMBM — Worked example
#
# =============================================================================
 
 
# -----------------------------------------------------------------------------
# LIBRARIES
# -----------------------------------------------------------------------------
library(sf)
library(terra)
library(ggplot2)
library(inlabru)
library(INLA)
library(sabinaNSDM)
library(sabinaJMBM)   # devtools::load_all(...)
 
 
# -----------------------------------------------------------------------------
# DATA
# -----------------------------------------------------------------------------
 
# Species name
SpeciesName <- "Quercus.petraea"
 
# Occurrence records
data(Quercus.petraea.xy.global, package = "sabinaNSDM")
data(Quercus.petraea.xy.regional, package = "sabinaNSDM")
 
# Environmental rasters (current)
data(expl.var.global, package = "sabinaNSDM")
data(expl.var.regional, package = "sabinaNSDM")
expl.var.global <- terra::unwrap(expl.var.global)
expl.var.regional <- terra::unwrap(expl.var.regional)
 
# Environmental rasters (new scenario)
data(new.env, package = "sabinaNSDM")
new.env <- terra::unwrap(new.env)
 
 
# -----------------------------------------------------------------------------
# DATA PREPARATION  (sabinaNSDM pipeline)
# -----------------------------------------------------------------------------
 
# Input object
myInput <- sabinaNSDM::NSDM.InputData(
  SpeciesName        = SpeciesName,
  spp.data.global    = Quercus.petraea.xy.global,
  spp.data.regional  = Quercus.petraea.xy.regional,
  expl.var.global    = expl.var.global,
  expl.var.regional  = expl.var.regional,
  new.env            = list(new.env),
  new.env.names      = "scenario1")
 
# Formatting
myFormatting <- sabinaNSDM::NSDM.FormattingData(
  myInput,
  nPoints            = 1000,
  Min.Dist.Global    = "resolution",
  Min.Dist.Regional  = "resolution",
  Background.method  = "random",
  save.output        = FALSE)
 
# Covariate selection
mySelvars <- sabinaNSDM::NSDM.SelectCovariates(
  myFormatting,
  maxncov.Global    = 3,
  maxncov.Regional  = 3,
  corcut            = 0.7,
  algorithms        = "glm",
  save.output       = FALSE)
 
 
# -----------------------------------------------------------------------------
# MESH  (spatial domain for the SPDE random field)
# -----------------------------------------------------------------------------
 
myMesh <- create_mesh(
  nsdm_obj         = mySelvars,
  edge             = c(2, 10),  # max triangle size c(inner, outer), in CRS units (here degrees)
  offset           = c(1, 5),   # domain extension c(inner, outer)
  buffer           = 0,
  boundary.method  = "raster_mask",
  remove_holes     = FALSE,
  proj.new.env     = TRUE,
  plot             = TRUE)       # set FALSE to skip the mesh plot
 
 
# -----------------------------------------------------------------------------
# JMBM MODEL FITTING
# -----------------------------------------------------------------------------
# coupling.intercept = how the regional intercept relates to the global one ("unpooled", "ordered_hierarchical", "bayesian_feedback")
# coupling.predictors = how shared covariates relate across scales ("unpooled", "ordered_hierarchical", "scale_decomposed", "bayesian_feedback")
# background.weights = corrects intercept bias when using presence-background data ("none" no correction, "auto" weights = 1 for presences, A/n_bg for background)
# inla.int.strategy = INLA hyperparameter integration ("eb", "ccd")

## ex covariate.effects
# cve <- list(global= list(bio12="linear", bio4=list(model="rw2", u=0.5, alpha=0.01)),
#             regional = list(bio1="linear", bio12="drop"), 
#             default="linear")
 
myModel <- JMBM.Modelling(
  jmbm_obj            = mySelvars,
  family              = binomial(link = "logit"),
  spde.mesh           = myMesh,
  regional.pcprior.range  = c(1.5, 0.05),   # S_re range prior
  regional.pcprior.sigma  = c(1.5, 0.01),   # S_re variance prior
  shared.pcprior.range = c(15, 0.05),    # S_shared range prior
  shared.pcprior.sigma = c(0.3, 0.01),   # S_shared variance prior
  coupling.intercept  = "ordered_hierarchical",
  coupling.predictors = "ordered_hierarchical",
  covariate.effects   = NULL,           # NULL = all covariates linear // cve 
  background.weights  = NULL,
  proj.new.env        = TRUE,
  cv.folds            = 1,               # 1 = no cross-validation
  n.threads           = 2,
  inla.int.strategy   = "eb",           # "ccd" for publication results
  seed                = 123,
  save.output         = FALSE)
 

# -----------------------------------------------------------------------------
# MODEL SUMMARY
# -----------------------------------------------------------------------------
summary(myModel)
 

# -----------------------------------------------------------------------------
# PLOTS
# -----------------------------------------------------------------------------
 
# Current suitability (posterior mean)
x11()
p_current <- plot(myModel, which = "pred", layer = "mean")
pts_reg <- sf::st_as_sf(
  mySelvars$SpeciesData.XY.Regional,
  coords = c("x", "y"),
  crs    = terra::crs(terra::unwrap(myModel$current.projections$pred), proj = TRUE))
p_current + ggplot2::geom_sf(data = pts_reg, colour = "black", size = 1.5, alpha = 0.4)
 
# Uncertainty (posterior standard deviation)
plot(myModel, which = "pred", layer = "sd")
 
# Residual spatial fields (if SPDE was fitted)
plot(myModel, which = "pred_Sre", layer = "mean")

# Broad-scale spatial field (Sshared)
plot(myModel, which = "pred_Sshared", layer = "mean")
 
# Future / alternative scenario
plot(myModel, which = "scenario1", layer = "mean")
