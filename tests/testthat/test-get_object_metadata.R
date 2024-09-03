test_that("object metadata retrieved", {
    expect_type(get_object_metadata(small_example_dataset), "list")
})
