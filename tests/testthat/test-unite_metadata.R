test_that("metadata united", {
    
    expect_error(
        unite_metadata(chevreul_sce, "nFeature_gene"),
        NA
    )
})
