test_that("plotting works", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
   chevreul_sce <- object_preprocess(chevreul_sce)
    expect_error(
        plot_violin(chevreul_sce, plot_var = "Mutation_Status", features = "Gene_0001"),
        NA
    )
})
