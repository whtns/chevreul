test_that("Read count metrics calculated", {
   chevreul_sce <- small_example_dataset
    expect_error(object_calcn(chevreul_sce), NA)
})
