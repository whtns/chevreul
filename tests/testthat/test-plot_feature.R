test_that("plotting works for sce", {
    expect_error(
        plot_feature(small_example_dataset, embedding = "UMAP", features = "Gene_0001"),
        NA
    )
})
