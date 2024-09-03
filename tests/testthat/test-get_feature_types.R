test_that("Feature type retrived", {
   chevreul_sce <- small_example_dataset

    expect_type(get_feature_types(chevreul_sce), "character")
})
