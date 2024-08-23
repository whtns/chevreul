test_that("percentage mito added", {
    expect_error(
        add_percent_mito(chevreul_sce),
        NA
    )
})
