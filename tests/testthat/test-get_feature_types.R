test_that("Feature type retrived", {
    

    expect_type(get_feature_types(chevreul_sce), "character")
})
