test_that("plotting works", {
    
    expect_error(
        plot_violin(chevreul_sce, plot_var = "batch", features = "NRL"),
        NA
    )
})
