test_that("plotting works", {
    
    expect_error(
        plot_var(chevreul_sce, "batch"),
        NA
    )
})
