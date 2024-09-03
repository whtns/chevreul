test_that("Metadata updated", {
   chevreul_sce <- small_example_dataset
    new_meta <- data.frame(row.names = colnames(chevreul_sce))
    new_meta$example <- "example"
    obj <- propagate_spreadsheet_changes(new_meta, chevreul_sce)
    expect_named(colData(obj), colnames(new_meta))
})
