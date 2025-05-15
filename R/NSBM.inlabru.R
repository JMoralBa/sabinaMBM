#' @name NSBM.inlabru
#'
#' @title Nested species distribution modeling ...
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling
#' with spatial structure via INLA and the inlabru.
#'
#' @param nsbm_obj An object of class `nsdm.vinput`, result from `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param mesh An INLA mesh object created externally with `create_mesh()`.
#' @param output A character. Either `"intensity"` or `"probability"` (default).
#' @param prior.range Numeric vector length 2. Prior on spatial range (e.g., `c(5, 0.01)`).
#' @param prior.sigma Numeric vector length 2. Prior on marginal standard deviation (e.g., `c(1, 0.01)`).
#'
#' @return A list with elements: `fit`, `pred`, and `pred_sp`.
#' @export
NSBM.inlabru <- function(nsbm_obj, mesh, output = "probability",
                         prior.range = c(5, 0.01),
                         prior.sigma = c(1, 0.01)) {

  if (!(output %in% c("probability", "intensity"))) {
    stop("Invalid 'output'. Use either 'probability' or 'intensity'.")
  }

  # Data preparation
  pp_regional <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Regional, coords = c("x", "y"))
  pp_global   <- sf::st_as_sf(nsbm_obj$SpeciesData.XY.Global, coords = c("x", "y"))
  sp_covglo   <- terra::unwrap(nsbm_obj$IndVar.Global.Selected)
  sp_covreg   <- terra::unwrap(nsbm_obj$IndVar.Regional.Selected)

  crs <- sf::st_crs(sp_covglo)
  sf::st_crs(pp_global) <- crs
  sf::st_crs(pp_regional) <- crs
  
  pp_global <- sf::st_transform(pp_global, crs)
  pp_regional <- sf::st_transform(pp_regional, crs)

  pts_reg <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pts_reg) <- crs
  pts_reg <- sf::st_transform(pts_reg, crs) 

  geom_comb <- c(sf::st_geometry(pp_global), sf::st_geometry(pts_reg))
  aux <- sf::st_sf(geometry = geom_comb)
  sf::st_crs(aux) <- crs
  
  # Spatial domain definition
  bdy_global <- sf::st_convex_hull(sf::st_union(aux))
  bdy_regional <- sf::st_union(sf::st_make_valid(sf::st_as_sf(raster::rasterToPolygons(raster::raster(sp_covreg)))))
  sf::st_crs(bdy_global) <- crs
  sf::st_crs(bdy_regional) <- crs

  # Mesh and SPDE specification
  matern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = prior.range,
    prior.sigma = prior.sigma
  )

  # Model components
  cmp_cov <- fcov(nsbm_obj, "sp_covglo", "sp_covreg")
  cmp <- as.formula(
    paste0("~ IGlobal(1) + IRegional(1) + spatial(geometry, model = matern) + ", cmp_cov$cmp)
  )

  if (output == "probability") {
    pp_regional$presence <- 1
    pseudo_regional <- sf::st_as_sf(nsbm_obj$Background.XY.Regional, coords = c("x", "y"))
    pseudo_regional$presence <- 0
    sf::st_crs(pseudo_regional) <- crs
    pseudo_regional <- sf::st_transform(pseudo_regional, crs)
    pp_regional <- rbind(pp_regional, pseudo_regional)

    pp_global$presence <- 1
    pseudo_global <- sf::st_as_sf(nsbm_obj$Background.XY.Global, coords = c("x", "y"))
    sf::st_crs(pseudo_global) <- crs
    pseudo_global$presence <- 0
    pseudo_global <- sf::st_transform(pseudo_global, crs) 
    pp_global <- rbind(pp_global, pseudo_global)

    # Likelihoods
    lik_global <- inlabru::like(
      family = "binomial",
      formula = as.formula(paste0("presence ~ IGlobal + spatial + ", cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = bdy_global,
      domain = list(geometry = mesh)
    )

    lik_regional <- inlabru::like(
      family = "binomial",
      formula = as.formula(paste0("presence ~ IRegional + spatial + ", cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = bdy_regional,
      domain = list(geometry = mesh)
    )

    pred_formula <- as.formula(
      paste0("~ 1 / (1 + exp(-(IRegional + spatial + ", cmp_cov$like$fregional, ")))")
    )

  } else {
    lik_global <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal + spatial + ", cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = bdy_global,
      domain = list(geometry = mesh)
    )

    lik_regional <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional + spatial + ", cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = bdy_regional,
      domain = list(geometry = mesh)
    )

    pred_formula <- as.formula(
      paste0("~ exp(IRegional + spatial + ", cmp_cov$like$fregional, ")")
    )
  }

  # Model fitting
  fit1 <- inlabru::bru(cmp, lik_global, lik_regional)

  # Prediction grid
  pred.df <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pred.df) <- crs
  pred.df <- sf::st_transform(pred.df, crs)

  # Predictions
  pred1 <- predict(fit1, pred.df, pred_formula) # distribución esperada 
  pred1_sp <- predict(fit1, pred.df, ~ spatial) # patrones espaciales residuales

  return(list(fit = fit1, pred = pred1, pred_sp = pred1_sp))
}


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

