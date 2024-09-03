test_that("plotting works for sce", {
    expect_error(
        plot_feature(chevreul_sce, embedding = "UMAP", features = "NRL"),
        NA
    )
})
