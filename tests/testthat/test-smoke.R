test_that("package loads and basic config constructors return lists", {
  cfg <- make_hazard_cfg()
  expect_true(is.list(cfg))

  sst_cfg <- make_sst_cfg(enabled = FALSE)
  expect_true(is.list(sst_cfg))
})
