test_that("SCE split", {
    expect_contains(
        obj <- splitByCol(small_example_dataset, "batch"),
        list()
    )
})
