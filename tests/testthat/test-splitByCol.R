test_that("SCE split", {
   chevreul_sce <- small_example_dataset
    expect_contains(
        obj <- splitByCol(chevreul_sce, "batch"),
        list()
    )
})
