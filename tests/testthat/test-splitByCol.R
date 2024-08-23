test_that("SCE split", {
    
    expect_contains(
        obj <- splitByCol(chevreul_sce, "batch"),
        list()
    )
})
