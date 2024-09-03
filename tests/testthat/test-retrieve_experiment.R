test_that("Assay retrieved", {
   mainExpName(small_example_dataset) <- "gene"
   retrieve_data <- retrieve_experiment(small_example_dataset, experiment = "gene")
    expect_equal(mainExpName(retrieve_data), "gene")
})
