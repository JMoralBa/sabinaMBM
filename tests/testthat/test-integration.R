# =============================================================================
# test-integration.R
# Integration tests: full mesh → fit → output structure.
# These tests run INLA — they are SLOW (~2–10 min each).
# Skipped automatically on CI. Run locally with: devtools::test()
# =============================================================================
 
# ── Shared setup (runs once per file, not per test) ──────────────────────────
# Uses real Fagus data from sabinaNSDM.
 
setup_fagus <- local({
  done <- FALSE
  mySelvars <- NULL
 
  function() {
    if(!done) {
      skip_if_not_installed("sabinaNSDM")
      data(Fagus.sylvatica.xy.global,   package = "sabinaNSDM", envir = environment())
      data(Fagus.sylvatica.xy.regional, package = "sabinaNSDM", envir = environment())
      data(expl.var.global,   package = "sabinaNSDM", envir = environment())
      data(expl.var.regional, package = "sabinaNSDM", envir = environment())
      data(new.env, package = "sabinaNSDM", envir = environment())
 
      expl.var.global   <- terra::unwrap(expl.var.global)
      expl.var.regional <- terra::unwrap(expl.var.regional)[[-3]]
      new.env           <- terra::unwrap(new.env)[[-3]]
 
      inp <- sabinaNSDM::NSDM.InputData(
        SpeciesName         = "Fagus.sylvatica",
        spp.data.global     = Fagus.sylvatica.xy.global,
        spp.data.regional   = Fagus.sylvatica.xy.regional,
        expl.var.global     = expl.var.global,
        expl.var.regional   = expl.var.regional,
        new.env             = list(new.env),
        new.env.names       = "scenario1")
 
      fmt <- sabinaNSDM::NSDM.FormattingData(
        inp, nPoints = 1000, Min.Dist.Global = "resolution",
        Min.Dist.Regional = "resolution",
        Background.method = "random", save.output = FALSE)
 
      mySelvars <<- sabinaNSDM::NSDM.SelectCovariates(
        fmt, maxncov.Global = 2, maxncov.Regional = 2,
        corcut = 0.7, algorithms = "glm", save.output = FALSE)
 
      done <<- TRUE
    }
    mySelvars
  }
})
 
setup_mesh <- local({
  done  <- FALSE
  mesh  <- NULL
 
  function(selvars) {
    if(!done) {
      mesh <<- create_mesh(
        nsdm_obj        = selvars,
        edge            = c(1, 2),
        offset          = c(0.5, 2),
        boundary.method = "convex_hull",
        plot            = FALSE)
      done <<- TRUE
    }
    mesh
  }
})
 
 
# ── Helper: check output structure ───────────────────────────────────────────
expect_valid_nsbm <- function(fit) {
  expect_s3_class(fit, "nsbm.inlabru")
  expect_named(fit, c("Species.Name", "args", "Selected.Variables.Global",
                       "Selected.Variables.Regional", "current.projections",
                       "new.projections", "marginals", "scale_params",
                       "pit_values", "diagnostic_plots", "Summary"),
               ignore.order = TRUE)
  expect_s4_class(terra::unwrap(fit$current.projections$pred), "SpatRaster")
  expect_true("mean" %in% names(terra::unwrap(fit$current.projections$pred)))
}
 
 
# ── Test 1: coupling.intercept = "unpooled" + coupling.predictors = NULL ─────
test_that("fit: unpooled intercept, no predictor coupling", {
  skip_on_ci()
  sv   <- setup_fagus()
  mesh <- setup_mesh(sv)
 
  fit <- NSBM.pure(
    nsbm_obj            = sv,
    family              = binomial(link = "logit"),
    spde.mesh           = mesh,
    local.pcprior.range  = c(2, 0.95),
    local.pcprior.sigma  = c(1, 0.01),
    coupling.intercept  = "ordered_hierarchical",
    coupling.predictors = NULL,
    proj.new.env        = TRUE,
    cv.folds = 1, n.threads = 1, seed = 1)
 
  expect_valid_nsbm(fit)
  # spatial SPDE is active (local.pcprior.range provided) → pred_local should exist
  expect_s4_class(terra::unwrap(fit$current.projections$pred_local), "SpatRaster")
  # latent SPDE is NOT active (shared.pcprior.range = NULL) → pred_shared should be NULL
  expect_null(fit$current.projections$pred_shared)
})
 
 
# ── Test 2: coupling.intercept = "ordered_hierarchical" + same for predictors ─
test_that("fit: ordered_hierarchical intercept + predictors", {
  skip_on_ci()
  sv   <- setup_fagus()
  mesh <- setup_mesh(sv)
 
  fit <- NSBM.pure(
    nsbm_obj            = sv,
    family              = binomial(link = "logit"),
    spde.mesh           = mesh,
    local.pcprior.range  = c(2, 0.95),
    local.pcprior.sigma  = c(1, 0.01),
    coupling.intercept  = "ordered_hierarchical",
    coupling.predictors = "ordered_hierarchical",
    proj.new.env        = TRUE,
    cv.folds = 1, n.threads = 1, seed = 1)
 
  expect_valid_nsbm(fit)
  # Spatial field should be present
  expect_s4_class(
    terra::unwrap(fit$current.projections$pred_local), "SpatRaster")
})
 
 
# ── Test 3: coupling.predictors = "scale_decomposed" ─────────────────────────
test_that("fit: scale_decomposed predictor coupling", {
  skip_on_ci()
  sv   <- setup_fagus()
  mesh <- setup_mesh(sv)
 
  fit <- NSBM.pure(
    nsbm_obj            = sv,
    family              = binomial(link = "logit"),
    spde.mesh           = mesh,
    local.pcprior.range  = c(2, 0.95),
    local.pcprior.sigma  = c(1, 0.01),
    coupling.intercept  = "unpooled",
    coupling.predictors = "scale_decomposed",
    proj.new.env        = FALSE,
    cv.folds = 1, n.threads = 1, seed = 1)
 
  expect_valid_nsbm(fit)
  # coefficients should include _ls or _ss suffix in the Summary fixed effects table
  fixed_names <- fit$Summary[["Fixed effects"]]$coef
  expect_true(any(grepl("_ls|_ss", fixed_names)))
})
 
 
# ── Test 4: coupling.intercept = NULL (regional-only) ────────────────────────
test_that("fit: regional-only model (coupling.intercept = NULL)", {
  skip_on_ci()
  sv   <- setup_fagus()
  mesh <- setup_mesh(sv)
 
  fit <- NSBM.pure(
    nsbm_obj            = sv,
    family              = binomial(link = "logit"),
    spde.mesh           = mesh,
    local.pcprior.range  = c(2, 0.95),
    local.pcprior.sigma  = c(1, 0.01),
    coupling.intercept  = NULL,
    coupling.predictors = NULL,
    proj.new.env        = FALSE,
    cv.folds = 1, n.threads = 1, seed = 1)
 
  expect_valid_nsbm(fit)
})
 
 
# ── Test 5: latent global SPDE ────────────────────────────────────────────────
test_that("fit: spatial + latent SPDE both active", {
  skip_on_ci()
  sv   <- setup_fagus()
  mesh <- setup_mesh(sv)
 
  fit <- NSBM.pure(
    nsbm_obj             = sv,
    family               = binomial(link = "logit"),
    spde.mesh            = mesh,
    local.pcprior.range   = c(2,  0.95),
    local.pcprior.sigma   = c(1,  0.01),
    shared.pcprior.range = c(8,  0.95),  # > 3× spatial range
    shared.pcprior.sigma = c(1,  0.01),
    coupling.intercept   = "unpooled",
    coupling.predictors  = NULL,
    proj.new.env         = FALSE,
    cv.folds = 1, n.threads = 1, seed = 1)
 
  expect_valid_nsbm(fit)
  expect_s4_class(
    terra::unwrap(fit$current.projections$pred_shared), "SpatRaster")
})
 
 
# ── Test 6: plot() does not throw ─────────────────────────────────────────────
test_that("plot.nsbm.inlabru: all which= options run without error", {
  skip_on_ci()
  sv   <- setup_fagus()
  mesh <- setup_mesh(sv)
 
  fit <- NSBM.pure(
    nsbm_obj            = sv,
    family              = binomial(link = "logit"),
    spde.mesh           = mesh,
    local.pcprior.range  = c(2, 0.95),
    local.pcprior.sigma  = c(1, 0.01),
    coupling.intercept  = "unpooled",
    coupling.predictors = NULL,
    proj.new.env        = TRUE,
    cv.folds = 1, n.threads = 1, seed = 1)
 
  expect_no_error(plot(fit, which = "pred"))
  expect_no_error(plot(fit, which = "pred_local"))
  expect_no_error(plot(fit, which = "scenario1"))
 
  # pred_shared absent → should give informative error
  expect_error(plot(fit, which = "pred_shared"), "pred_shared")
})
 
 
# ── Test 7: summary() does not throw ─────────────────────────────────────────
test_that("summary.nsbm.inlabru runs without error", {
  skip_on_ci()
  sv   <- setup_fagus()
  mesh <- setup_mesh(sv)
 
  fit <- NSBM.pure(
    nsbm_obj            = sv,
    family              = binomial(link = "logit"),
    spde.mesh           = mesh,
    local.pcprior.range  = c(2, 0.95),
    local.pcprior.sigma  = c(1, 0.01),
    coupling.intercept  = "unpooled",
    coupling.predictors = NULL,
    proj.new.env        = FALSE,
    cv.folds = 1, n.threads = 1, seed = 1)
 
  expect_no_error(summary(fit))
  expect_named(fit$Summary,
               c("Metadata", "Model fit", "Hyperparameters",
                 "Intercepts", "Fixed effects",
                 "Predictive performance", "Calibration & coverage", "Diagnostics"),
               ignore.order = TRUE)
})