test_that("preprocessing works", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    expect_error(
        object_preprocess(chevreul_sce),
        NA
    )
})
