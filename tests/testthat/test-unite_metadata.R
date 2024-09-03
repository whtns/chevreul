test_that("metadata united", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    expect_error(
    	unite_metadata(chevreul_sce, "Mutation_Status"),
        NA
    )
})
