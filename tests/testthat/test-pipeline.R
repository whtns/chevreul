test_that("integration pipeline works", {
    batches <- panc8 %>%
        Seurat::SplitObject(split.by = "tech")

    integrated_object <- object_integration_pipeline(batches)

    expect_equal(names(integrated_object@assays), c(names(panc8@assays), "integrated"))
})

test_that("object pipeline works", {
    processed_object <- object_pipeline(panc8)

    expect_named(processed_object@reductions, c("pca", "tsne", "umap"))
})
