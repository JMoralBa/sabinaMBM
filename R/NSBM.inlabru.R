#' @name NSBM.inlabru
#'
#' @title Nested species distribution modelling by means of hierarchical Bayesian models
#'
#' @description Estimation of presence probability by means of hierarchical 
#' Bayesian models.
#'
#' @param nsdm_obj An object of class \class{nsdm.vinput} with all the required data to fit the model.
#'
#' @param output A character. Either "intensity" or "probability".
#'
#' @return A list with three elements: \code{fit} (mode fit), \code{pred}
#' (prediction from the model) and \code{pred_sp} (prediction of the
#' spatial effect.
#'
#' @details TBC.
#'
#' @seealso \code{\link{NSDM.InputData}}, \code{\link{NSDM.FormattingData}}
#'
#' @examples
#' Example from GitHub
#' 
#' SpeciesName <- "Fagus.sylvativa"
#' 
#' # Species occurrences
#' data(Fagus.sylvatica.xy.global, package = "sabinaNSDM")
#' spp.data.global <- Fagus.sylvatica.xy.global
#' data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM")
#' spp.data.regional <- Fagus.sylvatica.xy.regional
#' 
#' data(expl.var.global, package = "sabinaNSDM")
#' data(expl.var.regional, package = "sabinaNSDM")
#' expl.var.global <- terra::unwrap(expl.var.global)
#' expl.var.regional <- terra::unwrap(expl.var.regional)
#'
#' # new escenarios
#' data(new.env, package = "sabinaNSDM")
#' new.env <- terra::unwrap(new.env)
#'
#' nsdm_input <- NSDM.InputData(SpeciesName = SpeciesName,
#'   spp.data.global = Fagus.sylvatica.xy.global,
#'   spp.data.regional = Fagus.sylvatica.xy.regional,
#'   expl.var.global = expl.var.global,
#'   expl.var.regional = expl.var.regional,
#'   new.env = new.env,
#'   new.env.names = "scenario1",
#'   Background.Global = NULL,
#'   Background.Regional = NULL,
#'   Absences.Global = NULL,
#'   Absences.Regional = NULL)
#' 
#' nsdm_finput <- NSDM.FormattingData(nsdm_input,
#'   nPoints = 100, # number of background points
#'   Min.Dist.Global = "resolution",
#'   Min.Dist.Regional = "resolution",
#'   Background.method = "random", # method “random" or "stratified” to generate background points 
#'   save.output = TRUE) #save outputs locally
#' 
#' nsdm_selvars <- NSDM.SelectCovariates(nsdm_finput,
#'   maxncov.Global = 3,   # Max number of covariates to be selected at the global scale
#'   maxncov.Regional = 3, # Max number of covariates to be selected at the regional scale
#'   corcut = 0.7, #  correlation threshold
#'   algorithms = c("glm"),
#'   ClimaticVariablesBands = NULL, # covariate bands to be excluded in the covariate selection at the regional scale
#'   save.output = TRUE)
#' 
#' # Probability of presence
#' nsdm_inlabru_prob <- NSBM.inlabru(nsdm_selvars)
#' # Intensity
#' nsdm_inlabru_int <- NSBM.inlabru(nsdm_selvars, "intensity")
#' 
#' @export
NSBM.inlabru <- function(nsbm_obj, output = "probability") {

  if (!(output %in% c("probability", "intensity"))) {
    stop("Wrong value for 'output'.")
  }

  # Define point patterns and covariates
  pp_regional <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Regional, coords = c("x", "y"))
  pp_global   <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Global, coords = c("x", "y"))
  sp_covglo   <- terra::unwrap(nsbm_obj$IndVar.Global.Selected)
  sp_covreg   <- terra::unwrap(nsbm_obj$IndVar.Regional.Selected)

  sf::st_crs(pp_global) <- sf::st_crs(sp_covglo)
  sf::st_crs(pp_regional) <- sf::st_crs(sp_covglo)

  # Define spatial boundaries
  aux <- sf::st_as_sf(sf::st_transform(sf::st_as_sfc(sf::st_as_sf(as.points(sp_covreg))), sf::st_crs(pp_global)))
  sf::st_geometry(aux) <- "geometry"
  aux <- rbind(pp_global, aux)
  bdy_global <- sf::st_convex_hull(sf::st_union(aux))
  bdy_regional <- sf::st_union(sf::st_make_valid(sf::st_as_sf(raster::rasterToPolygons(raster::raster(sp_covreg)))))
  bdy_global <- sf::st_buffer(bdy_global, 0.01)
  sf::st_crs(bdy_global) <- sf::st_crs(sp_covglo)
  sf::st_crs(bdy_regional) <- sf::st_crs(sp_covreg)

  # Define INLA mesh and SPDE
  mesh <- build_mesh(bdy_global, sf::st_crs(sp_covglo))
  matern <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(5, 0.01), prior.sigma = c(1, 0.01))

  # Define model components and formula
  cmp_cov <- fcov(nsbm_obj, "sp_covglo", "sp_covreg")
  cmp <- as.formula(paste0("~ IGlobal(1) + IRegional(1) + spatial(geometry, model = matern) + ", cmp_cov$cmp))

  # Define likelihoods and prediction formula
  if (output == "probability") {
    pp_regional$presence <- 1
    pseudo_regional <- sf::st_as_sf(nsbm_obj$Background.XY.Regional, coords = c("x", "y"))
    pseudo_regional$presence <- 0
    sf::st_crs(pseudo_regional) <- sf::st_crs(pp_regional)
    pp_regional <- rbind(pp_regional, pseudo_regional)

    pp_global$presence <- 1
    pseudo_global <- sf::st_as_sf(nsbm_obj$Background.XY.Global, coords = c("x", "y"))
    pseudo_global$presence <- 0
    sf::st_crs(pseudo_global) <- sf::st_crs(pp_global)
    pp_global <- rbind(pp_global, pseudo_global)

    lik_global <- inlabru::like(family = "binomial", formula = as.formula(paste0("presence ~ IGlobal + spatial + ", cmp_cov$like$fglobal)), data = pp_global, samplers = bdy_global, domain = list(geometry = mesh))
    lik_regional <- inlabru::like(family = "binomial", formula = as.formula(paste0("presence ~ IRegional + spatial + ", cmp_cov$like$fregional)), data = pp_regional, samplers = bdy_regional, domain = list(geometry = mesh))
    pred_formula <- as.formula(paste0("~ 1 / (1 + exp(-(IRegional + spatial + ", cmp_cov$like$fregional, ")))"))
  } else {
    lik_global <- inlabru::like(family = "cp", formula = as.formula(paste0("geometry ~ IGlobal + spatial + ", cmp_cov$like$fglobal)), data = pp_global, samplers = bdy_global, domain = list(geometry = mesh))
    lik_regional <- inlabru::like(family = "cp", formula = as.formula(paste0("geometry ~ IRegional + spatial + ", cmp_cov$like$fregional)), data = pp_regional, samplers = bdy_regional, domain = list(geometry = mesh))
    pred_formula <- as.formula(paste0("~ exp(IRegional + spatial + ", cmp_cov$like$fregional, ")"))
  }

  # Fit model and predict
  fit1 <- inlabru::bru(cmp, lik_global, lik_regional)
  pred.df <- sf::st_as_sf(as.points(sp_covreg))
  pred1 <- predict(fit1, pred.df, pred_formula)
  pred1_sp <- predict(fit1, pred.df, ~ spatial)

  return(list(fit = fit1, pred = pred1, pred_sp = pred1_sp))
}



  # Create bit of covariates for model formula
  # obj: Object with required data
  # spobjglo: Name of spatial object with layers (global)
  # spobjreg: Name of spatial object with layers (regional)
  fcov <- function(obj, spobjglo, spobjreg) {
    vars <- unique(c(obj$Selected.Variables.Global, obj$Selected.Variables.Regional))
    cmp1 <- paste(paste(vars, "(1)", sep = ""), collapse = " + ")

    cmpglobal <- paste(sapply(obj$Selected.Variables.Global, function(X) {
      paste0(X, "GL(main = ", spobjglo, ", main_layer = \"", X, "\", model = \"const\")")
    }), collapse = " + ")

    cmpregional <- paste(sapply(obj$Selected.Variables.Regional, function(X) {
      paste0(X, "RE(main = ", spobjreg, ", main_layer = \"", X, "\", model = \"const\")")
    }), collapse = " + ")

    fglobal <- paste(sapply(obj$Selected.Variables.Global, function(X) {
      paste0(X, " * ", X, "GL")
    }), collapse = " + ")

    fregional <- paste(sapply(obj$Selected.Variables.Regional, function(X) {
      paste0(X, " * ", X, "RE")
    }), collapse = " + ")

    list(
      cmp = paste(c(cmp1, cmpglobal, cmpregional), collapse = " + "),
      like = list(fglobal = fglobal, fregional = fregional)
    )
}


# Create INLA mesh for spatial modeling
# boundary: Outer boundary of the study area (sf object)
# crs: Coordinate reference system to assign to the mesh
# edge: Max edge lengths for mesh triangles (vector of 2 values)
build_mesh <- function(boundary, crs, edge = c(0.5, 1)) {
  mesh <- fmesher::fm_mesh_2d(
    boundary = boundary,
    max.edge = 4 * edge,
    offset = c(0.25, 0.5)
  )
  fmesher::fm_crs(mesh) <- crs
  return(mesh)
}
