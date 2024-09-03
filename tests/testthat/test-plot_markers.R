test_that("plotting works for sce", {
    expect_error(
        plot_markers(chevreul_sce, group_by = "gene_snn_res.0.2"),
        NA
    )
})
