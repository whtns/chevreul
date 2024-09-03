test_that("SCE split", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    expect_contains(
        obj <- splitByCol(chevreul_sce, "batch"),
        list()
    )
})
