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
  library(inlabru)
  library(sf)
  library(ggplot2)

  # Extraer elementos necesarios
  mesh <- nsbm_obj$mesh
  bdy_global <- nsbm_obj$bdy_global
  bdy_regional <- nsbm_obj$bdy_regional
  pp_global <- nsbm_obj$pp_global
  pp_regional <- nsbm_obj$pp_regional
  sp_covglo <- nsbm_obj$sp_covglo
  sp_covreg <- nsbm_obj$sp_covreg

  # Fórmulas de componentes y términos de covariables
  cmp_cov <- sabinaNSDM::fcov(nsbm_obj, "sp_covglo", "sp_covreg")

  cmp <- as.formula(paste0(
    "~ IGlobal(1) + IRegional(1) + spatial(geometry, model = matern) + ",
    cmp_cov$cmp
  ))

  # Preparar likelihoods y fórmula de predicción
  if (output == "probability") {
    pp_regional$presence <- 1
    pseudo_regional <- st_as_sf(nsbm_obj$Background.XY.Regional, coords = c("x", "y"))
    pseudo_regional$presence <- 0
    st_crs(pseudo_regional) <- st_crs(pp_regional)
    pp_regional <- rbind(pp_regional, pseudo_regional)

    pp_global$presence <- 1
    pseudo_global <- st_as_sf(nsbm_obj$Background.XY.Global, coords = c("x", "y"))
    pseudo_global$presence <- 0
    st_crs(pseudo_global) <- st_crs(pp_global)
    pp_global <- rbind(pp_global, pseudo_global)

    lik_global <- like(
      family = "binomial",
      formula = as.formula(paste0("presence ~ IGlobal + spatial + ", cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = st_as_sf(bdy_global),
      domain = list(geometry = mesh)
    )

    lik_regional <- like(
      family = "binomial",
      formula = as.formula(paste0("presence ~ IRegional + spatial + ", cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = st_as_sf(bdy_regional),
      domain = list(geometry = mesh)
    )

    pred_formula <- as.formula(paste0("~ 1 / (1 + exp(-(IRegional + spatial + ", cmp_cov$like$fregional, ")))"))

  } else {
    lik_global <- like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IGlobal + spatial + ", cmp_cov$like$fglobal)),
      data = pp_global,
      samplers = st_as_sf(bdy_global),
      domain = list(geometry = mesh)
    )

    lik_regional <- like(
      family = "cp",
      formula = as.formula(paste0("geometry ~ IRegional + spatial + ", cmp_cov$like$fregional)),
      data = pp_regional,
      samplers = st_as_sf(bdy_regional),
      domain = list(geometry = mesh)
    )

    pred_formula <- as.formula(paste0("~ exp(IRegional + spatial + ", cmp_cov$like$fregional, ")"))
  }

  # Ajustar modelo
  fit1 <- bru(cmp, lik_global, lik_regional)

  # Predicción espacial
  pred.df <- st_as_sf(as.points(sp_covreg))
  pred1 <- predict(fit1, pred.df, pred_formula)
  pred1_sp <- predict(fit1, pred.df, ~ spatial)

  return(list(
    fit = fit1,
    pred = pred1,
    pred_sp = pred1_sp
  ))
}

