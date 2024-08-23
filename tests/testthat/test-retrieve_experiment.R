test_that("Assay retrieved", {
    
    retrieve_data <- retrieve_experiment(chevreul_sce, experiment = "gene")
    expect_equal(mainExpName(retrieve_data), "gene")
})
