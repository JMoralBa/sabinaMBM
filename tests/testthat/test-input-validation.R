# =============================================================================
# test-input-validation.R
# Unit tests for all stop() / warning() checks in NSBM.pure()
# No INLA fit is triggered — these tests only reach the validation block.
# =============================================================================

test_that("nsbm_obj must be class nsdm.vinput", {
  expect_error(
    NSBM.pure(nsbm_obj = list(), family = binomial()),
    "nsdm.vinput"
  )
})


# ── family / link ─────────────────────────────────────────────────────────────

test_that("unsupported family throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj, family = Gamma(link = "inverse")),
    "Unsupported family"
  )
})

test_that("invalid link for binomial throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj, family = binomial(link = "log")),
    "link"
  )
})

test_that("family = 'cp' is accepted without error at validation stage", {
  obj  <- make_mock_vinput()
  mesh <- make_fake_mesh()
  # Should NOT throw a family/link error (may fail later for other reasons)
  expect_no_error_matching <- function(expr, pattern) {
    tryCatch(expr, error = function(e) {
      if(grepl(pattern, conditionMessage(e))) stop(conditionMessage(e))
    })
  }
  expect_no_error_matching(
    NSBM.pure(obj, family = "cp",
              spde.mesh = mesh,
              spde.pcprior.range = c(2, 0.95),
              spde.pcprior.sigma = c(1, 0.01)),
    "family|link"
  )
})


# ── coupling.intercept ────────────────────────────────────────────────────────

test_that("invalid coupling.intercept string throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj, coupling.intercept = "wrong_option"),
    "coupling.intercept"
  )
})

test_that("coupling.intercept = NULL with no regional covariates throws error", {
  obj <- make_mock_vinput(vg = c("bio1"), vr = character(0))
  expect_error(
    NSBM.pure(obj, coupling.intercept = NULL),
    "regional"
  )
})

test_that("coupling.intercept != NULL with no global covariates throws error", {
  obj <- make_mock_vinput_regional_only()
  mesh <- make_fake_mesh()
  expect_error(
    NSBM.pure(obj,
              coupling.intercept = "unpooled",
              spde.mesh = mesh,
              spde.pcprior.range = c(2, 0.95),
              spde.pcprior.sigma = c(1, 0.01)),
    "global component"
  )
})

test_that("coupling.intercept = NULL with latent SPDE throws error", {
  obj  <- make_mock_vinput(vg = character(0), vr = c("bio1"))
  mesh <- make_fake_mesh()
  expect_error(
    NSBM.pure(obj,
              coupling.intercept   = NULL,
              spde.mesh            = mesh,
              latent.pcprior.range = c(5, 0.95),
              latent.pcprior.sigma = c(1, 0.01)),
    "Latent global SPDE"
  )
})


# ── SPDE / mesh ───────────────────────────────────────────────────────────────

test_that("SPDE priors with NULL mesh throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj,
              spde.mesh          = NULL,
              spde.pcprior.range = c(2, 0.95),
              spde.pcprior.sigma = c(1, 0.01),
              coupling.intercept = NULL),  # evita el check de mesh previo
    "spde.mesh"
  )
})

test_that("mesh of wrong class throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj,
              spde.mesh          = list(),   # not inla.mesh
              spde.pcprior.range = c(2, 0.95),
              spde.pcprior.sigma = c(1, 0.01)),
    "inla.mesh"
  )
})

test_that("latent range < 3x spatial range triggers warning", {
  # Requiere mesh real — cubierto en test-integration.R
  # "fit: spatial + latent SPDE both active"
  skip("Covered in integration tests")
})


# ── covariate.effects ─────────────────────────────────────────────────────────

test_that("covariate.effects must be a list", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj, covariate.effects = "linear"),
    "must be a list"
  )
})

test_that("covariate.effects with invalid top-level key throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj, covariate.effects = list(wrong_key = "linear")),
    "Invalid entries"
  )
})

test_that("covariate.effects default = 'rw2' as string throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj, covariate.effects = list(default = "rw2")),
    "must be 'linear' or 'drop'"
  )
})

test_that("rw2 without u/alpha throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj,
              covariate.effects = list(
                global = list(bio1 = list(model = "rw2", u = 0.5))  # missing alpha
              )),
    "alpha"
  )
})

test_that("global RW2 with active latent SPDE throws error", {
  obj  <- make_mock_vinput()
  mesh <- make_fake_mesh()
  expect_error(
    NSBM.pure(obj,
              spde.mesh            = mesh,
              spde.pcprior.range   = c(2,  0.95),
              spde.pcprior.sigma   = c(1,  0.01),
              latent.pcprior.range = c(10, 0.95),
              latent.pcprior.sigma = c(1,  0.01),
              covariate.effects    = list(
                global = list(bio1 = list(model = "rw2", u = 0.5, alpha = 0.01))
              )),
    "RW2 global effects"
  )
})


# ── coupling.predictors ───────────────────────────────────────────────────────

test_that("invalid coupling.predictors string throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj, coupling.predictors = "wrong_mode"),
    "coupling.predictors"
  )
})

test_that("coupling.predictors list with invalid key throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj,
              coupling.predictors = list(wrong_key = "unpooled")),
    "Invalid entries"
  )
})

test_that("coupling.predictors list with invalid variable mode throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj,
              coupling.predictors = list(
                default   = "unpooled",
                variables = list(bio1 = "bad_mode")
              )),
    "Invalid coupling mode"
  )
})

test_that("ordered_hierarchical on variable missing in global throws error", {
  # bio4 is regional-only
  obj <- make_mock_vinput(vg = c("bio1"), vr = c("bio1", "bio4"))
  expect_error(
    NSBM.pure(obj,
              coupling.predictors = list(
                default   = "unpooled",
                variables = list(bio4 = "ordered_hierarchical")  # bio4 not in global
              )),
    "Missing in global"
  )
})

test_that("mixing bayesian_feedback with ordered_hierarchical throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj,
              coupling.intercept  = "bayesian_feedback",
              coupling.predictors = "ordered_hierarchical"),
    "Cannot mix"
  )
})


# ── seed / n.threads ──────────────────────────────────────────────────────────

test_that("non-numeric seed throws error", {
  obj <- make_mock_vinput()
  expect_error(
    NSBM.pure(obj, seed = "abc"),
    "seed"
  )
})
