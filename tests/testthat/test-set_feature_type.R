test_that("Feature type set", {
    new_object <- set_feature_type(small_example_dataset, "gene")
    expect_equal(mainExpName(new_object), "gene")
})
