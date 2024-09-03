test_that("returns a character vector", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    expect_type(
        metadata_from_object(chevreul_sce),
        "character"
    )
})
