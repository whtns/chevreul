test_that("percentage mito added", {
	chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    expect_error(
        add_percent_mito(chevreul_sce),
        NA
    )
})
