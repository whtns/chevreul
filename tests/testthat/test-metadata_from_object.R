test_that("returns a character vector", {
   chevreul_sce <- small_example_dataset
    expect_type(
        metadata_from_object(chevreul_sce),
        "character"
    )
})
