test_that("percentage mito added", {
	chevreul_sce <- small_example_dataset
    expect_error(
        add_percent_mito(chevreul_sce),
        NA
    )
})
