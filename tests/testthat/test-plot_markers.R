test_that("plotting works for sce", {
    expect_error(
        plot_markers(small_example_dataset, group_by = "gene_snn_res.1"),
        NA
    )
})
