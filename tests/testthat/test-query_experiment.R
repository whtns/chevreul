test_that("possible to query experiment names", {
   mainExpName(small_example_dataset) <- "gene"
    expect_equal(query_experiment(small_example_dataset, "gene"), TRUE)
})
