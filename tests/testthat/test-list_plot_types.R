test_that("List of variables to plot produced", {
    var_list <- list_plot_types(small_example_dataset)
    expect_contains(var_list, list())
})
