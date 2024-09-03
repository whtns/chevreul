test_that("Feature type retrived", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)

    expect_type(get_feature_types(chevreul_sce), "character")
})
