# 12/05/2025

#dir
 setwd("C:/Users/Jennifer Morales/Documents/UAM-INLABRU/")

# R version 4.3.2 (2023-10-31 ucrt)

# Required packages
 library(sabinaNSDM)
 library(glmnet)
 library(stringi)
 library(sf)
 library(terra)
 #library("raster")
 library(covsel)
 library(biomod2)
 library(ecospat)
 library(fs)
 library(sgsR)
 library(remotes)
 #remotes::install_github("inlabru-org/inlabru", dependencies = TRUE)
 library("inlabru")
 library("INLA")
 library("tidyterra")
 library("ggplot2")
 library(rnaturalearth)
 library(rnaturalearthdata)
 library(ggtext)
 library(concaveman)
 #library(mgcv)
 library(pROC)
 library(covsel)

 devtools::load_all("C:/Users/Jennifer Morales/Documents/UAM-INLABRU/repo/sabinaINLA") # path al paquete clonado en pc


## EXAMPLE INLABRU
## Data
 SpeciesName <- "Fagus.sylvativa"
 
# Species occurrences
 data(Fagus.sylvatica.xy.global, package = "sabinaNSDM")
 spp.data.global <- Fagus.sylvatica.xy.global
 data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM")
 spp.data.regional <- Fagus.sylvatica.xy.regional
 
 data(expl.var.global, package = "sabinaNSDM")
 data(expl.var.regional, package = "sabinaNSDM")
 expl.var.global <- terra::unwrap(expl.var.global)
 expl.var.regional <- terra::unwrap(expl.var.regional)
 expl.var.regional <- expl.var.regional[[-3]] #rm sand

 # new escenarios
 data(new.env, package = "sabinaNSDM")
 new.env <- terra::unwrap(new.env)
 new.env <- new.env[[-3]] #rm sand
 #new.env2 <- r_con_isla # sale del script prueba_different_env.R


## INPUT DATA
 myInput <- sabinaNSDM::NSDM.InputData(
   SpeciesName = SpeciesName,
   spp.data.global = Fagus.sylvatica.xy.global,
   spp.data.regional = Fagus.sylvatica.xy.regional,
   expl.var.global = expl.var.global,
   expl.var.regional = expl.var.regional,
   new.env = new.env,
   new.env.names = "scenario1",
   Background.Global = NULL,
   Background.Regional = NULL,
   Absences.Global = NULL,
   Absences.Regional = NULL)

## FORMATING DATA
 myFormatting <- sabinaNSDM::NSDM.FormattingData(
   myInput,
   nPoints = 300,
   Min.Dist.Global = "resolution",
   Min.Dist.Regional = "resolution",
   Background.method = "random", 
   save.output = FALSE)

## SELECT VARIABLES
 mySelvars <- sabinaNSDM::NSDM.SelectCovariates(
   myFormatting,
   maxncov.Global = 3,
   maxncov.Regional = 3,
   corcut = 0.7,d
   algorithms = c("glm"),
   ClimaticVariablesBands = NULL,
   save.output = FALSE)

 #save(mySelvars, file = "repo/sabinaINLA/data/mySelvars.RData")
 #save(mySelvars, file = "repo/sabinaINLA/data/mySelvars_diff_env.RData")
 #load("repo/sabinaINLA/data/mySelvars.RData")

## sabinaINLA
## MESH
# check crs para ser coiherente con los valores de edge, offset y buffer (grados, metros?).
 crs(terra::unwrap(mySelvars$IndVar.Global.Selected))

 myMesh <- create_mesh(
   nsdm_obj = mySelvars,  # Objeto nsdm.vinput de sabunaNSDM
   edge = c(0.5, 1), #c(1, 2),        # tamaño máx de triángulos c(interior, exterior)
   offset = c(0.25, 0.5), #c(2, 3),      # expansión del dominio c(interior, exterior) 
   buffer = 0,            # margen adicional del borde  
   boundary.method = "convex_hull", # "convex_hull", "raster_mask", "cancave_hull"
   concavity = NULL,      # Numeric que ajusta inner mesh a los puntos en concave_hull
   remove_holes = FALSE,  # filtra los holes de inner mesh para raster_mask
   proj.new.env = FALSE,  # extiende la mesh a new.env si necesario cuando raster_mask
   plot = TRUE)           # Plotea la malla
 

## INLA PURE HIERARCHICAL
 myPred.pure <- NSBM.pure(
   nsbm_obj = mySelvars,       # Objeto nsdm.vinput de sabinaNSDM
   family = binomial(link="logit"),
   spde.mesh = myMesh,              # Add efecto esapcial SPDE wirh mesh of create_mesh()
   spde.pcprior.range = c(5, 0.95),   # Prior para el rango espacial: c(valor, prob. de ser menor)
   spde.pcprior.sigma = c(1, 0.01),   # Prior para la desviación estándar: c(valor, prob. de ser mayor)
   latent.pcprior.range = c(5, 0.95), 
   latent.pcprior.sigma = c(1, 0.01),
   nested.intercept = TRUE,    # TRUE: IRegional se modela como desviación de IGlobal
   covariate.pcprior.smoothness = NULL,  #NULL, list(u=0.5, alpha=0.01) or "auto". Not running!!!
   proj.new.env = FALSE,       # Project new scenarios
   cv.folds = 1,               # k-folds para cross-validation
   n.threads =1,               # hilos de inla/inlabru
   seed = NULL,
   save.output = FALSE) 
 
 summary(myPred.pure)
 #str(myPred.pure)
 #save(myPred.pure, file = "myPred_pure.RData")
 #rm(myPred.pure)
 #load("myPred_pure.RData")
 #?NSBM.pure
 #names((myPred.pure))

 # plot
 pred <- terra::unwrap(myPred.pure$current.projections$pred)
 #pred_sp <- terra::unwrap(myPred.pure$current.projections$pred_sp)
 #new.proj <- terra::unwrap(myPred.pure$new.projections[[1]])

 #names(pred)
 #terra::plot(pred[[1]])

 x11()
 map.prob <- ggplot() +
   geom_raster(data = pred, aes(x = x, y = y, fill = mean)) +
   scale_fill_distiller(palette = "Spectral", name = "Suitability", na.value = "transparent") +
   labs(x = "Longitude", y = "Latitude", title = "Pure hierarchical") +
   theme_minimal()
 points_pres <- sf::st_as_sf(mySelvars$SpeciesData.XY.Regional, coords = c("x", "y"))
 sf::st_crs(points_pres) <- sf::st_crs(pred)
 map.prob <- map.prob +
   geom_sf(data = points_pres, color = "black", size = 1.5, alpha = 0.4)
 map.prob


## INLA MULTIPLY
 myPred.multiply <- NSBM.multiply(
   nsbm_obj = mySelvars,
   output = "probability",
   family = "binomial",
   link = "logit",
   spde.mesh = myMesh,
   spde.pcprior.range = c(5, 0.01),
   spde.pcprior.sigma = c(1, 0.01),
   method = "geometric",
   rescale = TRUE,
   proj.new.env = TRUE,
   cv.folds = 1, 
   n.threads = 1,
   seed = NULL,
   save.output = FALSE,
   save.independent = FALSE)

 str(myPred.multiply)
 summary(myPred.multiply)
 #?NSBM.multiply

 pred <- terra::unwrap(myPred.multiply$current.projections$pred.multiply)

 x11()
 map.multiply <- ggplot() +
   geom_raster(data = pred, aes(x = x, y = y, fill = Fagus.sylvativa.Current)) +
   scale_fill_distiller(palette = "Spectral", name = "Suitability", na.value = "transparent") +
   labs(x = "Longitude", y = "Latitude", title = "Multiply") +
   theme_minimal()
 # Añadir puntos de presencia
 points_pres <- sf::st_as_sf(mySelvars$SpeciesData.XY.Regional, coords = c("x", "y"))
 sf::st_crs(points_pres) <- sf::st_crs(pred)
 map.multiply <- map.multiply +
   geom_sf(data = points_pres, color = "black", size = 1.5, alpha = 0.4)
 map.multiply


## INLA COVARIATE
 myPred.covariate <- NSBM.covariate(
   nsbm_obj = mySelvars,
   family = binomial(link="logit"),
   spde.mesh = myMesh,
   spde.pcprior.range = c(5, 0.01),
   spde.pcprior.sigma = c(1, 0.01),
   rm.corr = TRUE,
   corcut = 0.7,         # correlation threshold
   proj.new.env = TRUE,
   cv.folds = 1, 
   n.threads = 1,
   seed = NULL,
   save.output = TRUE,
   save.independent = FALSE)

 str(myPred.covariate)
 summary(myPred.covariate)
 #?NSBM.covariate

 pred <- terra::unwrap(myPred.covariate$current.projections$pred.covariate)

 x11()
 map.covariate <- ggplot() +
   geom_raster(data = pred, aes(x = x, y = y, fill = mean)) +
   scale_fill_distiller(palette = "Spectral", name = "Suitability", na.value = "transparent") +
   labs(x = "Longitude", y = "Latitude", title = "Covariate") +
   theme_minimal()
 points_pres <- sf::st_as_sf(mySelvars$SpeciesData.XY.Regional, coords = c("x", "y"))
 sf::st_crs(points_pres) <- sf::st_crs(pred)
 map.covariate <- map.covariate +
   geom_sf(data = points_pres, color = "black", size = 1.5, alpha = 0.4)
 map.covariate

