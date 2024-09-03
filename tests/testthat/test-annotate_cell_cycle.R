test_that("Annotation of readcount works", {
    data(cc.genes.cyclone)
    expect_error(
        annotate_cell_cycle(small_example_dataset),
        NA
    )
})
