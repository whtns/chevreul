context("cluster methods")

test_that("seurat preprocessed", {
  processed_seurat <- seurat_preprocess(panc8$gene)

  # variable features calculated
  expect_gt(length(Seurat::VariableFeatures(processed_seurat)), 0)

  # data scaled
  expect_false(all(is.na(slot(processed_seurat$RNA, "scale.data"))))
})

test_that("clustering workflow works", {
  clustered_seu <- clustering_workflow(panc8)

})
