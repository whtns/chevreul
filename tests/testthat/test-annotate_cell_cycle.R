test_that("Annotation of readcount works", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    data(cc.genes.cyclone)
    expect_error(
        annotate_cell_cycle(chevreul_sce),
        NA
    )
})
