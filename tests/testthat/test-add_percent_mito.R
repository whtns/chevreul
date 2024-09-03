test_that("percentage mito added", {
    expect_error(
        add_percent_mito(small_example_dataset),
        NA
    )
})
