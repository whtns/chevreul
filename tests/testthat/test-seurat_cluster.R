context("cluster methods")

test_that("object preprocessed", {
  processed_object <- object_preprocess(panc8$gene)

  # variable features calculated
  expect_gt(length(Seurat::VariableFeatures(processed_object)), 0)

  # data scaled
  expect_false(all(is.na(slot(processed_object[["gene"]], "scale.data"))))
})

test_that("clustering workflow works", {
  clustered_object <- clustering_workflow(panc8)
})
