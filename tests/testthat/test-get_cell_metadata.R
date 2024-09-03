test_that("Metadata pulled", {
   chevreul_sce <- small_example_dataset
    meta <- get_cell_metadata(chevreul_sce)
    expect_contains(colnames(meta), colnames(colData(chevreul_sce)))
})
