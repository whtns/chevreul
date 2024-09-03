test_that("returns a character vector", {
    expect_type(
        metadata_from_object(small_example_dataset),
        "character"
    )
})
