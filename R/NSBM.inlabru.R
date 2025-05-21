#' @name NSBM.inlabru
#'
#' @title Nested species distribution modeling ...
#'
#' @description Fits a nested species distribution model (NSDM) using hierarchical Bayesian modeling
#' with spatial structure via INLA and inlabru.
#'
#' @param nsbm_obj An object of class `nsdm.vinput`, result from `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param mesh An INLA mesh object created externally with `create_mesh()`. Required only if \code{spatial = TRUE}. Ignored if \code{spatial = FALSE}.
#' @param output A character `"intensity"` or `"probability"` (default).
#' @param prior.range Numeric vector length 2. Prior on spatial range (e.g., `c(5, 0.01)`).
#' @param prior.sigma Numeric vector length 2. Prior on marginal standard deviation (e.g., `c(1, 0.01)`).
#' @param spatial Logical. Include spatial latent field (SPDE) in the model (default: TRUE). 
#' @param proj.new.env Logical. Whether to compute predictions under new scenarios (default: TRUE).
#' @param seed Optional integer. If provided, sets a random seed for reproducibility.
#'
#' @return A named list of class `nsbm.inlabru` with the following elements:
#' \item{Species.Name}{Species name}
#' \item{args}{List of arguments used in the model fitting.}
#' \item{Selected.Variables.Global}{Names of selected global-scale covariates.}
#' \item{Selected.Variables.Regional}{Names of selected regional-scale covariates.}
#' \item{current.projections}{List with: fitted model (`fit`), prediction (`pred`), spatial field (`pred_sp`).}
#' \item{new.projections}{List of projections to new.env (if `proj.new.env = TRUE`).}
#' \item{Summary}{Empty placeholder (data.frame) for future model summaries.}
#'
#' @export
NSBM.inlabru <- function(nsbm_obj, 
                         mesh = NULL, 
                         output = "probability",
                         spatial = TRUE,
                         prior.range = c(5, 0.01),
                         prior.sigma = c(1, 0.01),
                         proj.new.env = TRUE,
                         seed = NULL) {

  if(!inherits(nsbm_obj, "nsdm.vinput")) {
    stop("The 'nsbm_obj' must be of class 'nsdm.vinput'.")
  }
  if(!(output %in% c("probability", "intensity"))) {
    stop("Invalid 'output'. Use either 'probability' or 'intensity'.")
  }
  if(spatial && is.null(mesh)) {
    stop("If spatial = TRUE, you must provide a mesh object using create_mesh().")
  }
  if(!is.null(seed)){
    set_seed(seed)
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

  # Description SPDE
  if(spatial) {
    matern <- INLA::inla.spde2.pcmatern(
      mesh,
      prior.range = prior.range,
      prior.sigma = prior.sigma
    )
  }

  # Model components
  cmp_cov <- fcov(nsbm_obj, "sp_covglo", "sp_covreg")
  cmp_formula <- if(spatial) {
    paste0("~ IGlobal(1) + IRegional(1) + spatial(geometry, model = matern) + ", cmp_cov$cmp)
  } else {
    paste0("~ IGlobal(1) + IRegional(1) + ", cmp_cov$cmp)
  }

  cmp <- as.formula(cmp_formula)

  if(output == "probability") {
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
    f_spatial <- if(spatial) " + spatial" else ""

    lik_global <- inlabru::like(
      family = "binomial",
      formula = as.formula(paste0("presence ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = bdy_global,
      domain = list(geometry = mesh)
    )

    lik_regional <- inlabru::like(
      family = "binomial",
      formula = as.formula(paste0("presence ~ IRegional", f_spatial, " + ", cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = bdy_regional,
      domain = list(geometry = mesh)
    )

    pred_formula <- if(spatial) {
      as.formula(paste0("~ 1 / (1 + exp(-(IRegional + spatial + ", cmp_cov$like$fregional, ")))"))
      } else {
        as.formula(paste0("~ 1 / (1 + exp(-(IRegional + ", cmp_cov$like$fregional, ")))"))
      }

  } else {
    lik_global <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal", f_spatial, " + ", cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = bdy_global,
      domain = list(geometry = mesh)
    )

    lik_regional <- inlabru::like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional", f_spatial, " + ", cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = bdy_regional,
      domain = list(geometry = mesh)
    )

    pred_formula <- if(spatial) {
      as.formula(paste0("~ exp(IRegional + spatial + ", cmp_cov$like$fregional, ")"))
    } else {
      as.formula(paste0("~ exp(IRegional + ", cmp_cov$like$fregional, ")"))
    }
  }

  # Model fitting
  fit <- inlabru::bru(cmp, lik_global, lik_regional)

  pred.df <- sf::st_as_sf(as.points(sp_covreg))
  sf::st_crs(pred.df) <- crs
  pred.df <- sf::st_transform(pred.df, crs)

  # Predictions
  pred <- predict(fit, pred.df, pred_formula) # distribución esperada 
  pred_sp <- if(spatial) {  # patrones espaciales residuales
    predict(fit, pred.df, ~ spatial)
  } else {
    NULL
  } 

  # Project new scenarios
  proj_list <- list()
  if (proj.new.env && !is.null(nsbm_obj$Scenarios)) {
    for (sc in names(nsbm_obj$Scenarios)) {
      scen_rast <- terra::unwrap(nsbm_obj$Scenarios[[sc]])
      scen_df <- sf::st_as_sf(as.points(scen_rast))
      sf::st_crs(scen_df) <- crs
      scen_df <- sf::st_transform(scen_df, crs)
      proj_pred <- predict(fit, scen_df, pred_formula)
      proj_list[[paste0("proj_", sc)]] <- proj_pred
    }
  }

  sabina <- list(
      Species.Name = nsbm_obj$SpeciesName,
      args = list(
        output = output,
        spatial = spatial,
        prior.range = prior.range,
        prior.sigma = prior.sigma
      ),
      Selected.Variables.Global = nsbm_obj$Selected.Variables.Global,
      Selected.Variables.Regional = nsbm_obj$Selected.Variables.Regional,
      current.projections = list(
        fit = fit,
        pred = pred,
        pred_sp = pred_sp
      ),
      new.projections = proj_list,
      Summary = data.frame() #@@@JMB pendiente
  )
 
  attr(sabina, "class") <- "nsbm.inlabru"
  return(sabina)

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


