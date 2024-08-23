test_that("Dimensional reduction successful", {
    
    reduce_obj <- object_reduce_dimensions(chevreul_sce)
    expect_contains(reducedDimNames(reduce_obj), c("PCA", "TSNE", "UMAP"))
})
