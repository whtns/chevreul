test_that("metadata united", {
    expect_error(
        unite_metadata(small_example_dataset, "Mutation_Status"),
        NA
    )
})
