test_that("Metadata pulled", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    meta <- get_cell_metadata(chevreul_sce)
    expect_contains(colnames(meta), colnames(colData(chevreul_sce)))
})
