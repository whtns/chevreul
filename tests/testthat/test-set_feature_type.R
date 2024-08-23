test_that("Feature type set", {
    
    new_object <- set_feature_type(chevreul_sce, "transcript")
    expect_equal(mainExpName(new_object), "transcript")
})
