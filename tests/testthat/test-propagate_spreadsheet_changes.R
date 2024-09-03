test_that("Metadata updated", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    new_meta <- data.frame(row.names = colnames(chevreul_sce))
    new_meta$example <- "example"
    obj <- propagate_spreadsheet_changes(new_meta, chevreul_sce)
    expect_named(colData(obj), colnames(new_meta))
})
