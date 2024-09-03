test_that("Metadata pulled", {
    meta <- get_cell_metadata(small_example_dataset)
    expect_contains(colnames(meta), colnames(colData(small_example_dataset)))
})
