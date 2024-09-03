test_that("Assay retrieved", {
   chevreul_sce <- small_example_dataset
   mainExpName(chevreul_sce) <- "gene"
   retrieve_data <- retrieve_experiment(chevreul_sce, experiment = "gene")
    expect_equal(mainExpName(retrieve_data), "gene")
})
