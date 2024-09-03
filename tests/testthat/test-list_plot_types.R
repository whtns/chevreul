test_that("List of variables to plot produced", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    var_list <- list_plot_types(chevreul_sce)
    expect_contains(var_list, list())
})
