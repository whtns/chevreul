test_that("possible to query experiment names", {
   chevreul_sce <- small_example_dataset
   mainExpName(chevreul_sce) <- "gene"
    expect_equal(query_experiment(chevreul_sce, "gene"), TRUE)
})
