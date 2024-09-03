test_that("List of variables to plot produced", {
   chevreul_sce <- small_example_dataset
    var_list <- list_plot_types(chevreul_sce)
    expect_contains(var_list, list())
})
