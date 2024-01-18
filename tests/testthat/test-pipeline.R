test_that("integration pipeline works", {
    batches <- panc8 %>%
        Seurat::SplitObject(split.by = "tech")

    integrated_seu <- seurat_integration_pipeline(batches)

    expect_equal(names(integrated_seu@assays), c(names(panc8@assays), "integrated"))
})

test_that("seurat pipeline works", {
    processed_seu <- seurat_pipeline(panc8)

    expect_named(processed_seu@reductions, c("pca", "tsne", "umap"))
})
