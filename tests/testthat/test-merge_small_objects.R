test_that("SCE objects merged", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    expect_contains(obj <- merge_small_objects(
        "small_batch1" = chevreul_sce[, seq(100)],
        "small_batch2" = chevreul_sce[, seq(101, 200)]
    ), list())
})
