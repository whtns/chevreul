test_that("Read count metrics calculated", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    expect_error(object_calcn(chevreul_sce), NA)
})
