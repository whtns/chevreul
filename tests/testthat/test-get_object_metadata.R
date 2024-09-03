test_that("object metadata retrieved", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    expect_type(get_object_metadata(chevreul_sce), "list")
})
