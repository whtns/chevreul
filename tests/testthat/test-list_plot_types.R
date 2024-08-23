test_that("List of variables to plot produced", {
    
    var_list <- list_plot_types(chevreul_sce)
    expect_contains(var_list, list())
})
