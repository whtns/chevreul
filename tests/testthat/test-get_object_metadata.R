test_that("object metadata retrieved", {
   chevreul_sce <- small_example_dataset
    expect_type(get_object_metadata(chevreul_sce), "list")
})
