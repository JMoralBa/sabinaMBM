# =============================================================================
# test-resolve-helpers.R
# Unit tests for .resolve_covariate_effects() and .resolve_coupling_predictor()
# These are pure functions — no rasters, no INLA required.
# =============================================================================

# ── .resolve_covariate_effects ────────────────────────────────────────────────

test_that("NULL covariate.effects returns linear by default", {
  r <- sabinaNSBM:::.resolve_covariate_effects("bio1", "global", NULL)
  expect_equal(r$model, "linear")
})

test_that("explicit 'linear' string returns linear", {
  ce <- list(global = list(bio1 = "linear"))
  r  <- sabinaNSBM:::.resolve_covariate_effects("bio1", "global", ce)
  expect_equal(r$model, "linear")
})

test_that("explicit 'drop' string returns drop", {
  ce <- list(global = list(bio1 = "drop"))
  r  <- sabinaNSBM:::.resolve_covariate_effects("bio1", "global", ce)
  expect_equal(r$model, "drop")
})

test_that("valid rw2 list returns model=rw2 with u and alpha", {
  ce <- list(global = list(bio1 = list(model = "rw2", u = 0.5, alpha = 0.01)))
  r  <- sabinaNSBM:::.resolve_covariate_effects("bio1", "global", ce)
  expect_equal(r$model, "rw2")
  expect_equal(r$u,     0.5)
  expect_equal(r$alpha, 0.01)
})

test_that("rw2 as bare string throws error", {
  ce <- list(global = list(bio1 = "rw2"))
  expect_error(
    sabinaNSBM:::.resolve_covariate_effects("bio1", "global", ce),
    "RW2 must be expressed as a list"
  )
})

test_that("rw2 list without alpha throws error", {
  ce <- list(global = list(bio1 = list(model = "rw2", u = 0.5)))
  expect_error(
    sabinaNSBM:::.resolve_covariate_effects("bio1", "global", ce),
    "alpha"
  )
})

test_that("invalid keyword in covariate.effects throws error", {
  ce <- list(global = list(bio1 = "smooth"))
  expect_error(
    sabinaNSBM:::.resolve_covariate_effects("bio1", "global", ce),
    "Invalid keyword"
  )
})

test_that("variable not listed inherits from default", {
  ce <- list(default = "drop")
  r  <- sabinaNSBM:::.resolve_covariate_effects("bio99", "global", ce)
  expect_equal(r$model, "drop")
})

test_that("variable not listed, no default, returns linear", {
  ce <- list(global = list(bio1 = "linear"))
  r  <- sabinaNSBM:::.resolve_covariate_effects("bio99", "global", ce)
  expect_equal(r$model, "linear")
})

test_that("scale-specific entry overrides default", {
  ce <- list(
    global  = list(bio1 = "drop"),
    default = "linear"
  )
  r <- sabinaNSBM:::.resolve_covariate_effects("bio1", "global", ce)
  expect_equal(r$model, "drop")
})


# ── .resolve_coupling_predictor ───────────────────────────────────────────────

test_that("NULL coupling.predictors returns 'unpooled'", {
  r <- sabinaNSBM:::.resolve_coupling_predictor("bio1", NULL, vg = c("bio1"))
  expect_equal(r, "unpooled")
})

test_that("character mode is returned for variable in global", {
  r <- sabinaNSBM:::.resolve_coupling_predictor(
    "bio1", "ordered_hierarchical", vg = c("bio1"))
  expect_equal(r, "ordered_hierarchical")
})

test_that("hierarchical mode silently falls back to unpooled for regional-only var", {
  # bio4 is NOT in vg → should silently downgrade to unpooled
  r <- sabinaNSBM:::.resolve_coupling_predictor(
    "bio4", "ordered_hierarchical", vg = c("bio1"))
  expect_equal(r, "unpooled")
})

test_that("list form: variable-specific mode overrides default", {
  cp <- list(default = "unpooled", variables = list(bio1 = "scale_decomposed"))
  r  <- sabinaNSBM:::.resolve_coupling_predictor("bio1", cp, vg = c("bio1"))
  expect_equal(r, "scale_decomposed")
})

test_that("list form: unlisted variable uses default", {
  cp <- list(default = "ordered_hierarchical", variables = list(bio12 = "unpooled"))
  r  <- sabinaNSBM:::.resolve_coupling_predictor("bio1", cp, vg = c("bio1"))
  expect_equal(r, "ordered_hierarchical")
})

test_that("explicit hierarchical on regional-only var is NOT silently downgraded", {
  # When explicitly set, it should NOT downgrade — the check in NSBM.pure catches this
  cp <- list(default = "unpooled", variables = list(bio4 = "ordered_hierarchical"))
  r  <- sabinaNSBM:::.resolve_coupling_predictor("bio4", cp, vg = c("bio1"))
  expect_equal(r, "ordered_hierarchical")  # explicit → stays; NSBM.pure will stop()
})
