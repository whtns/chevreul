test_that("Metadata updated", {
    new_meta <- data.frame(row.names = colnames(small_example_dataset))
    new_meta$example <- "example"
    obj <- propagate_spreadsheet_changes(new_meta, small_example_dataset)
    expect_named(colData(obj), colnames(new_meta))
})
