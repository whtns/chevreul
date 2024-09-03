test_that("Annotation of readcount works", {
   chevreul_sce <- small_example_dataset
    data(cc.genes.cyclone)
    expect_error(
        annotate_cell_cycle(chevreul_sce),
        NA
    )
})
