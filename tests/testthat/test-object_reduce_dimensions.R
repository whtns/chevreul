test_that("Dimensional reduction successful", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
   chevreul_sce <- object_preprocess(chevreul_sce)
   
    reduce_obj <- object_reduce_dimensions(chevreul_sce)
    expect_contains(reducedDimNames(reduce_obj), c("PCA", "TSNE", "UMAP"))
})
