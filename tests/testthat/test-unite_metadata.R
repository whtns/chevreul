test_that("metadata united", {
   chevreul_sce <- small_example_dataset
    expect_error(
    	unite_metadata(chevreul_sce, "Mutation_Status"),
        NA
    )
})
