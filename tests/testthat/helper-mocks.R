# =============================================================================
# helper-mocks.R
# Minimal mock objects for unit tests — no sabinaNSDM or INLA required.
# Loaded automatically by testthat before every test file.
# =============================================================================

# Dependencies needed by integration tests (sabinaNSDM pipeline)
library(glmnet)
library(stringi)
library(covsel)
library(sabinaNSDM)

# ── Tiny synthetic rasters (10×10 pixels, WGS84) ─────────────────────────────
.make_raster <- function(varnames, nrow = 10, ncol = 10) {
  r <- terra::rast(
    nrows = nrow, ncols = ncol,
    xmin = -5, xmax = 5, ymin = 35, ymax = 45,
    crs  = "EPSG:4326"
  )
  stk <- lapply(varnames, function(v) {
    vals     <- r
    terra::values(vals) <- stats::rnorm(nrow * ncol)
    names(vals) <- v
    vals
  })
  do.call(c, stk)
}

# ── Tiny presence-background data.frame ──────────────────────────────────────
.make_xy <- function(n = 20) {
  data.frame(
    x    = stats::runif(n, -4, 4),
    y    = stats::runif(n, 36, 44),
    resp = sample(c(0L, 1L), n, replace = TRUE)
  )
}

# ── Minimal nsdm.vinput mock ──────────────────────────────────────────────────
#' Build a minimal nsdm.vinput object.
#' @param vg  Character vector: global covariate names
#' @param vr  Character vector: regional covariate names
#' @param with_scenarios Logical: include Scenarios slot
make_mock_vinput <- function(vg   = c("bio1", "bio12"),
                             vr   = c("bio1", "bio4"),
                             with_scenarios = FALSE) {

  sp_glo <- .make_raster(vg)
  sp_reg <- .make_raster(vr)

  obj <- list(
    Species.Name                  = "Mock.species",
    Selected.Variables.Global     = vg,
    Selected.Variables.Regional   = vr,
    IndVar.Global.Selected        = if(length(vg) > 0) terra::wrap(sp_glo) else NULL,
    IndVar.Regional.Selected      = if(length(vr) > 0) terra::wrap(sp_reg) else NULL,
    SpeciesData.XY.Global         = .make_xy(30),
    SpeciesData.XY.Regional       = .make_xy(20),
    Response.Global               = NULL,
    Response.Regional             = NULL,
    Scenarios                     = if(with_scenarios) list(scenario1 = terra::wrap(sp_reg)) else NULL
  )
  class(obj) <- "nsdm.vinput"
  obj
}

# ── Convenience: mock with NO global covariates ───────────────────────────────
make_mock_vinput_regional_only <- function() {
  make_mock_vinput(vg = character(0), vr = c("bio1", "bio4"))
}

# ── Fake inla.mesh (minimal S3 object to pass class check) ───────────────────
make_fake_mesh <- function() {
  m <- list()
  class(m) <- "inla.mesh"
  m
}
