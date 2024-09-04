test_that("Feature names obained", {
    names <- get_features(small_example_dataset)
    expect_equal(names, rownames(small_example_dataset))
})
