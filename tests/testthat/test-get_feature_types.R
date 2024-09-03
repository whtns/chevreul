test_that("Feature type retrived", {

    expect_type(get_feature_types(small_example_dataset), "character")
})
