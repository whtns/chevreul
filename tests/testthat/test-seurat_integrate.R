test_that("integration doesn't drift", {
  integrated_panc8 <- reintegrate_object(panc8)

  reintegrated_panc8 <- reintegrate_object(integrated_panc8)

  expect_equal(integrated_panc8@reductions$umap, reintegrated_panc8@reductions$umap)
})
